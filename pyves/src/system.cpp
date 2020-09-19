#include "system.hpp"


namespace _pyves
{
    void System::setThreads(std::size_t num)
    {
        // threads = std::max(num, static_cast<std::size_t>(2));
        threads = num;

#ifdef PYVES_USE_TBB
	    tbb::global_control control(tbb::global_control::max_allowed_parallelism, threads);

        if(threads == 1)
        {
            task_arena.initialize(threads, 1);
        }
        else
        {
            task_arena.initialize(threads);
        }

#else
        executor.reset(new tf::Executor(threads));
#endif
    }



    bool System::particleIsFree(const Particle& subject) const
    {
        return std::all_of(std::begin(particles), std::end(particles), [&](const auto& p)
        {
            if(std::isnan(subject.sigma) || std::isnan(p.sigma))
            {
                throw std::logic_error("System::particleIsFree found particle where sigma is NaN");
            }
            return subject == p ? true : box.distance(subject.position, p.position) > 1.1224f*(subject.sigma+p.sigma)/2;
        });
    }



    bool System::particleIsFree(const Particle& subject, REAL cutoff) const
    {
        const REAL cutoff_squared = cutoff*cutoff;
        return std::all_of(std::begin(particles), std::end(particles), [&](const auto& p)
        { 
            return subject == p ? true : box.squaredDistance(subject.position, p.position) > cutoff_squared;
        });
    }



    std::size_t System::numParticlesInCells() const
    {
        return std::accumulate(std::begin(cells), std::end(cells), std::size_t(0), [](auto i, const Cell& c) { return i + c.particles.size(); });
    }



    bool System::assertIntegrity()
    {
        return all(
            numParticlesInCells() == particles.size(),
            std::all_of(std::begin(particles), std::end(particles), [](const auto& p) { return p.assertIntegrity(); }),
            std::all_of(std::begin(cells), std::end(cells), [](auto& c) { return c.assertIntegrity(); }),
            std::isfinite(threads),
            std::isfinite(temperature)
        );
    }



    void System::singleSimulationStep()
    {
#ifdef PYVES_USE_TBB
        task_arena.execute([&]
        {
            prepareSimulationStep();
            applyToCells([&](const Cell& c){ cellStep(c);});
        });
#else
        prepareSimulationStep();
        applyToCells([&](const Cell& c){ cellStep(c);});
#endif
        ++_cell_update_step_count;
        ++_neighbor_update_step_count;
    }



    void System::multipleSimulationSteps(const unsigned long _max)
    {
        for(unsigned long step = 0; step < _max; ++step )
        {
            singleSimulationStep();
        }
    }
    


    void System::shuffle()
    {
        // std::cout << __func__ << " " << internal_step_count << "  for tranlsation " << translation_alignment() << "\n";
#ifdef PYVES_USE_TBB
        tbb::parallel_for_each(std::begin(cells), std::end(cells), [&] (Cell& cell) 
#else
        tf::Taskflow taskflow;
        taskflow.for_each_static(std::begin(cells), std::end(cells), [&] (Cell& cell) 
#endif
        {
            cell.shuffle(); 
        });
#ifndef PYVES_USE_TBB
        executor->run(taskflow).get();
#endif
    }



    void System::reorderCells()
    {
        // std::cout << __func__ << " " << internal_step_count << "  for tranlsation " << translation_alignment() << "\n";
#ifdef PYVES_USE_TBB
        tbb::parallel_for_each(std::begin(cells), std::end(cells), [&] (Cell& cell) 
#else
        tf::Taskflow taskflow;
        taskflow.for_each_static(std::begin(cells), std::end(cells), [&] (Cell& cell) 
#endif
        {
            // std::cout << cell.repr() << "\n";
            auto leavers = cell.particlesOutOfBounds();
            
            for(Particle& leaver : leavers)
            {
                bool was_added = false;
                // std::cout << cell.proximity.size() << "\n";
                for(Cell& proximity_cell : cell.proximity)
                {
                    // std::cout << " " << proximity_cell.repr() << "\n";
                    was_added = proximity_cell.try_add(leaver);
                    if(was_added) 
                    {
                        break;
                    }
                }
                if(!was_added)
                {
                    throw std::logic_error("Particle out of bound was not added to another cell\n" + leaver.repr() +" and " + cell.repr());
                }
                assert(was_added);
                cell.removeParticle(leaver);
                assert(!cell.contains(leaver));
            }
            // cell.state = (cell.particles.size() > 0) ? CellState::IDLE : CellState::FINISHED;
            // cell.shuffle(); 
        });
#ifndef PYVES_USE_TBB
        executor->run(taskflow).get();
#endif
    }



    void System::makeNeighborLists()
    {
        // std::cout << __func__ << " " << internal_step_count << "  for tranlsation " << translation_alignment() << "\n";
#ifdef PYVES_USE_TBB
        tbb::parallel_for_each(std::begin(cells), std::end(cells), [&] (Cell& cell) 
#else
        tf::Taskflow taskflow;
        taskflow.for_each(std::begin(cells), std::end(cells), [&] (Cell& cell) 
#endif
        {
            cell.updateRegionParticles();
            for(Particle& p : cell.particles)
            {
                p.updateNeighborList(cell.region_particles, box, interaction_cutoff + 1);
            }
        });
#ifndef PYVES_USE_TBB
        executor->run(taskflow).get();
#endif
    }



    void System::prepareSimulationStep()
    {
        if( _cell_update_step_count >= cell_update_interval )
        {
            reorderCells();
            shuffle();
            _cell_update_step_count = 0;
        }
        if( _neighbor_update_step_count >= neighbor_update_interval )
        {
            makeNeighborLists();
            _neighbor_update_step_count = 0;
        }
    }


    
    REAL System::potentialEnergy() const
    {
        REAL sum = 0;
        for(std::size_t i = 0; i < particles.size(); ++i)
        {
            sum += particles[i].potentialEnergy(box, interaction_cutoff);
        }
        return sum/2;
    }


    
    REAL System::potentialEnergyConcurrent()
    {
#ifdef PYVES_USE_TBB
        REAL sum = 0;
        task_arena.execute([&]{
            using range_type = tbb::blocked_range<ParticleContainer::iterator>;
            sum = tbb::parallel_reduce(
                range_type(std::begin(particles), std::end(particles)), 
                REAL(0), 
                [&](const range_type& r, REAL init ){ 
                    return std::accumulate( r.begin(), r.end(), init, [&](REAL _v, const Particle& p){
                        return _v + p.potentialEnergy(box, interaction_cutoff);
                    }); 
                }, 
                std::plus<REAL>()
            );
        });
        return sum/2;
#else 
        REAL sum = 0;
        tf::Taskflow taskflow;
        taskflow.transform_reduce_static(
            std::begin(particles), 
            std::end(particles), 
            sum, 
            std::plus<REAL>(), 
            [&](const Particle& p){
                return p.potentialEnergy(box, interaction_cutoff); 
            }
        );
        executor->run(taskflow).get();
        return sum/2;
#endif
    }



    void System::benchmark(std::size_t num)
    {
        auto timeFuncInvocation = [](std::size_t num, auto&& func, auto&&... params) {
            // get time before function invocation
            const auto& start = std::chrono::high_resolution_clock::now();
            // function invocation using perfect forwarding
            for(volatile std::size_t i = 0; i<num; ++i)
                std::forward<decltype(func)>(func)(std::forward<decltype(params)>(params)...);
            // get time after function invocation
            const auto& stop = std::chrono::high_resolution_clock::now();
            return std::chrono::duration<double>(stop - start).count();
        };

        std::cout << "num x potentialEnergy                " <<  timeFuncInvocation(num, [&]{potentialEnergy();})<< " s\n"; 
        std::cout << "num x potentialEnergyConcurrent      " <<  timeFuncInvocation(num, [&]{potentialEnergyConcurrent();}) << " s\n";
        std::cout << "num x singleSimulationStep           " <<  timeFuncInvocation(num, [&]{singleSimulationStep();}) << " s\n";
        std::cout << "multipleSimulationSteps(num)         " <<  timeFuncInvocation(1, [&]{multipleSimulationSteps(num);}) << " s\n";
    }



    void System::cellStep(const Cell& cell)
    {
        // TODO:FIXME:TODO:FIXME: CAUTION: this function is exactly what it should look like.
        // do not change order or modify any function call

        // std::cout << cell.repr() << "\n";
        REAL last_energy_value;
        REAL energy_after;
        std::uniform_real_distribution<REAL> dist_coords(-translation_alignment(),translation_alignment());
        std::uniform_real_distribution<REAL> dist_orientation(-rotation_alignment(),rotation_alignment());
        CARTESIAN translation;//(dist_coords(pseudo_engine), dist_coords(pseudo_engine), dist_coords(pseudo_engine));
        CARTESIAN random_vector;//(dist_coords(pseudo_engine), dist_coords(pseudo_engine), dist_coords(pseudo_engine));
        CARTESIAN orientation_before;//(dist_coords(pseudo_engine), dist_coords(pseudo_engine), dist_coords(pseudo_engine));

        std::vector<ParticleRefContainer::const_iterator> iterators(cell.particles.size());
        std::iota(std::begin(iterators), std::end(iterators), std::cbegin(cell.particles));
        std::shuffle(std::begin(iterators), std::end(iterators), RandomEngine.pseudo_engine);

        for(const auto& iterator : iterators)
        {
            Particle& particle = *iterator;
            
            // coordinates move
            {
                do
                {
                    translation = CARTESIAN::Random() * dist_coords(RandomEngine.pseudo_engine);
                }
                while(translation.squaredNorm() > translation_alignment()*translation_alignment());

                // last_energy_value = cell.potentialEnergy(particle, interaction_cutoff);
                last_energy_value = particle.potentialEnergy(box, interaction_cutoff);

                if(particle.trySetPosition(particle.position+translation))
                {
                    // energy_after = cell.potentialEnergy(particle, interaction_cutoff);
                    energy_after = particle.potentialEnergy(box, interaction_cutoff);                

                    // rejection
                    if(!acceptByMetropolis(energy_after - last_energy_value, temperature))
                    {
                        translation_alignment.rejected();
                        particle.position -= translation;
                        // particle.position = particle.position - translation;
                    }
                    // acctance
                    else
                    {
                        translation_alignment.accepted();
                        last_energy_value = energy_after;
                    }
                }
                else
                {
                    translation_alignment.rejected();
                }
            }

        // TODO:FIXME:TODO:FIXME: CAUTION: this function is exactly what it should look like.
        // do not change order or modify any function call

            do
            {
                random_vector = CARTESIAN::Random();
            }
            // to not oversample with points outside of the unity sphere
            while(random_vector.squaredNorm() > 1);

            orientation_before = particle.orientation; 
            // particle.setOrientation(Eigen::AngleAxis<REAL>(dist_orientation(RandomEngine.pseudo_engine), random_vector) * particle.getOrientation());

            if(particle.trySetOrientation(Eigen::AngleAxis<REAL>(dist_orientation(RandomEngine.pseudo_engine), random_vector) * particle.getOrientation()))
            {
                // energy_after = cell.potentialEnergy(particle, interaction_cutoff);
                energy_after = particle.potentialEnergy(box, interaction_cutoff);

                // rejection
                if(!acceptByMetropolis(energy_after - last_energy_value, temperature))
                {
                    rotation_alignment.rejected();
                    particle.setOrientation(orientation_before);
                }
                // acctance
                else
                {
                    rotation_alignment.accepted();
                    last_energy_value = energy_after;
                }
            }
            else
            {
                rotation_alignment.rejected();
            }
        }

        // TODO:FIXME:TODO:FIXME: CAUTION: this function is exactly what it should look like.
        // do not change order or modify any function call
    }



    void System::makeInteractionLookupTable(ParticleContainer unqiues)
    {
        lookup_table.clear();
        lookup_table = calculateInteractionLookupTable(unqiues);
    }



    void bind_system(py::module& m)
    {
	    py::bind_vector<ParticleContainer>( m, "ParticleContainer" )
            .def(py::init<>())
            // .def("clear", &ParticleContainer::clear)
            // .def("pop_back", &ParticleContainer::pop_back)
            .def("__len__", [](const ParticleContainer& v) { return v.size(); })
            .def("__iter__", [](ParticleContainer& v) 
            {
                return py::make_iterator(std::begin(v), std::end(v));
            }, py::keep_alive<0, 1>())
            .def("__repr__", [](const ParticleContainer& v) {
                return "ParticleContainer\n[\n"
                    + std::accumulate(std::begin(v), std::end(v), std::string(""), [](std::string s, const Particle& p) 
                    { 
                        return s+"\t"+p.repr()+",\n";
                    })
                    + "]";
            })
        ;



	    py::bind_vector<CellContainer>(m, "CellContainer" )
            .def(py::init<>())
            // .def("clear", &ParticleContainer::clear)
            // .def("pop_back", &ParticleContainer::pop_back)
            .def("__len__", [](const CellContainer& v) { return v.size(); })
            .def("__iter__", [](CellContainer& v) 
            {
                return py::make_iterator(std::begin(v), std::end(v));
            }, py::keep_alive<0, 1>())
            .def("__repr__", [](const CellContainer& v) {
                return "CellContainer\n[\n"
                    + std::accumulate(std::begin(v), std::end(v), std::string(""), [](std::string s, const Cell& p) 
                    { 
                        return s+"\t"+p.repr()+",\n";
                    })
                    + "]";
            })
        ;

        py::bind_map<LookupTable_t>(m, "LookupTable");


        py::class_<System>(m, "System", py::dynamic_attr())
            .def(py::init<>())
            .def_property("threads", [](const System& s){ return s.threads; }, &System::setThreads)
            .def_readwrite("temperature", &System::temperature)
            .def_readwrite("box", &System::box)
            .def_readwrite("particles", &System::particles)
            .def_readwrite("interaction_cutoff", &System::interaction_cutoff)
            .def_readwrite("cells", &System::cells)
            .def_readwrite("cell_update_interval", &System::cell_update_interval)
            .def_readwrite("neighbor_update_interval", &System::neighbor_update_interval)
            .def("particleIsFree", static_cast<bool (System::*)(const Particle&) const>(&System::particleIsFree))
            .def("particleIsFree", static_cast<bool (System::*)(const Particle&, REAL) const>(&System::particleIsFree))
            .def("shuffle", &System::shuffle)
            .def("numParticlesInCells", &System::numParticlesInCells)
            .def_readwrite("translation", &System::translation_alignment)
            .def_readwrite("rotation", &System::rotation_alignment)
            .def("assertIntegrity", &System::assertIntegrity)
            .def("prepareSimulationStep", &System::prepareSimulationStep)
            .def("singleSimulationStep", &System::singleSimulationStep)
            .def("multipleSimulationSteps", &System::multipleSimulationSteps)
            .def("makeLookupTableFrom", &System::makeInteractionLookupTable)
            .def("makeNeighborLists", &System::makeNeighborLists)
            .def("potentialEnergy", &System::potentialEnergy)
            .def("potentialEnergyConcurrent", &System::potentialEnergyConcurrent)
            .def("benchmark", &System::benchmark)
        ;
    }
}