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



    bool System::particleIsFree(const Particle& subject, REAL cutoff) const
    {
        const REAL cutoff_squared = cutoff*cutoff;
        return (interaction::external_potential(subject, box, interaction_surface_width, interaction_cutoff) < 1e-3) && 
            std::all_of(std::begin(particles), std::end(particles), [&](const Particle& p)
            { 
                return subject == p ? true : box.squaredDistance(subject.position, p.position) > cutoff_squared;
            }
        );
    }



    bool System::particleIsFree(const Particle& subject) const
    {
        // return particleIsFree(subject, 1.1224f*(subject.sigma+p.sigma)/2);
        return (interaction::external_potential(subject, box, interaction_surface_width, interaction_cutoff) < 1e-3) and 
            std::all_of(std::begin(particles), std::end(particles), [&](const Particle& p)
            {
                if(std::isnan(subject.sigma) || std::isnan(p.sigma))
                {
                    throw std::logic_error("System::particleIsFree found particle where sigma is NaN");
                }
                return subject == p ? true : box.distance(subject.position, p.position) > 1.1224f*(subject.sigma+p.sigma)/2;
            }
        );
    }



    std::size_t System::numParticlesInCells() const
    {
        return std::accumulate(std::begin(cells), std::end(cells), std::size_t(0), [](auto i, const Cell& c) { return i + c.particles.size(); });
    }



    bool System::assertIntegrity()
    {
        return all(
            numParticlesInCells() == particles.size(),
            std::all_of(std::begin(particles), std::end(particles), [](const Particle& p) { return p.assertIntegrity(); }),
            std::all_of(std::begin(cells), std::end(cells), [](Cell& c) { return c.assertIntegrity(); }),
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
        globalExchange();

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
        if(!executor)
        {
            setThreads(1);
        }
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



    // void System::reorderCell(Cell& cell)
    // {
    //     // std::cout << cell.repr() << "\n";
    //     auto leavers = cell.particlesOutOfBounds();
        
    //     for(Particle& leaver : leavers)
    //     {
    //         bool was_added = false;
    //         // std::cout << cell.proximity.size() << "\n";
    //         for(Cell& proximity_cell : cell.proximity)
    //         {
    //             // std::cout << " " << proximity_cell.repr() << "\n";
    //             was_added = proximity_cell.try_add(leaver);
    //             if(was_added) 
    //             {
    //                 break;
    //             }
    //         }
    //         if(!was_added)
    //         {
    //             throw std::logic_error("Particle out of bound was not added to another cell\n" + leaver.repr() +" and " + cell.repr());
    //         }
    //         assert(was_added);
    //         cell.removeParticle(leaver);
    //         assert(!cell.contains(leaver));
    //     }
    //     // cell.state = (cell.particles.size() > 0) ? CellState::IDLE : CellState::FINISHED;
    //     // cell.shuffle(); 
    // }



    void System::reorderCells()
    {
        // std::cout << __func__ << " " << internal_step_count << "  for tranlsation " << translation_alignment() << "\n";
#ifdef PYVES_USE_TBB
        tbb::parallel_for_each(std::begin(cells), std::end(cells), [&] (Cell& cell) 
#else
        if(!executor)
        {
            setThreads(1);
        }
        tf::Taskflow taskflow;
        taskflow.for_each_static(std::begin(cells), std::end(cells), [&] (Cell& cell) 
#endif
        {
            cell.reorder();
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
        if(!executor)
        {
            setThreads(1);
        }
        tf::Taskflow taskflow;
        taskflow.for_each(std::begin(cells), std::end(cells), [&] (Cell& cell) 
#endif
        {
            cell.updateRegionParticles();
            for(Particle& p : cell.particles)
            {
                p.updateNeighborList(cell.region_particles, box, interaction_cutoff + 0.1);
            }
        });
#ifndef PYVES_USE_TBB
        executor->run(taskflow).get();
#endif
    }



    void System::prepareSimulationStep()
    {
        if( _cell_update_step_count >= cell_update_interval)
        {
            reorderCells();
            shuffle();
            _cell_update_step_count = 0;
        }
        if( _neighbor_update_step_count >= neighbor_update_interval)
        {
            makeNeighborLists();
            _neighbor_update_step_count = 0;
        }
    }


    
    REAL System::totalEnergy(const Particle& p) const
    {
        return p.potentialEnergy(box, interaction_cutoff) + interaction::surface_potential(p, box, interaction_surface_width, interaction_cutoff) * (interaction_surface ? 1 : 0);
    }


    
    REAL System::totalEnergy()
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
                        return _v + totalEnergy(p);
                    }); 
                }, 
                std::plus<REAL>()
            );
        });
        return sum/2;
#else 
        if(!executor)
        {
            setThreads(1);
        }
        REAL sum = 0;
        tf::Taskflow taskflow;
        taskflow.transform_reduce_static(
            std::begin(particles), 
            std::end(particles), 
            sum, 
            std::plus<REAL>(), 
            [&](const Particle& p){
                return totalEnergy(p); 
            }
        );
        executor->run(taskflow).get();
        return sum/2;
#endif
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


    
    Eigen::Matrix<REAL,Eigen::Dynamic, 1> System::particleEnergies() const
    {
        Eigen::Matrix<REAL,Eigen::Dynamic, 1> energies(particles.size());
        for(std::size_t i = 0; i < particles.size(); ++i)
        {
            energies(i) = particles[i].potentialEnergy(box, interaction_cutoff)/2;
        }
        return energies;
    }


    
    Eigen::Matrix<REAL,Eigen::Dynamic, 1> System::particleChiValues() const
    {
        Eigen::Matrix<REAL,Eigen::Dynamic, 1> chis(particles.size());
        for(std::size_t i = 0; i < particles.size(); ++i)
        {
            chis(i) = particles[i].chi(box, interaction_cutoff)/2;
        }
        return chis;
    }


    
    Eigen::Matrix<REAL,Eigen::Dynamic, 1> System::particleSurfacePotentialValues() const
    {
        Eigen::Matrix<REAL,Eigen::Dynamic, 1> values(particles.size());
        for(std::size_t i = 0; i < particles.size(); ++i)
        {
            values(i) = interaction::surface_potential(particles.at(i), box, interaction_surface_width, interaction_cutoff);
            if(values(i) < 1e5)
            {
                auto z = particles.at(i).getz();
                while(z > box.getLengthZ())
                {
                    z -= box.getLengthZ();
                }
                if(
                    interaction_surface and 
                    z < interaction_cutoff and
                    particles.at(i).surface_affinity_translation > 1e-3
                )
                throw std::logic_error("impossible surface potential = " + std::to_string(values(i)) 
                                       + " for " + particles.at(i).repr() + " z=" +std::to_string(z)
                                       + " with surface_width " + std::to_string(interaction_surface_width)
                                       + " and cutoff " + std::to_string(interaction_cutoff) );
            }
        }
        return values;
    }


    
    Eigen::Matrix<REAL,Eigen::Dynamic, 1> System::particleExternalPotentialValues() const
    {
        Eigen::Matrix<REAL,Eigen::Dynamic, 1> values(particles.size());
        for(std::size_t i = 0; i < particles.size(); ++i)
        {
            values(i) = interaction::external_potential(particles.at(i), box, interaction_surface_width, interaction_cutoff);
        }
        return values;
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
        if(!executor)
        {
            setThreads(1);
        }
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

        std::cout << "\nBENCHMARK of energy functions with " << threads << " threads:\n";
        const auto single = timeFuncInvocation(num, [&]{potentialEnergy();});
        const auto concurrent = timeFuncInvocation(num, [&]{potentialEnergyConcurrent();});
        std::cout << num << " x potentialEnergy                " << single << " s\n"; 
        std::cout << num << " x potentialEnergyConcurrent      " << concurrent << " s\n";
        std::cout << "scaling factor: " << single/(concurrent*threads) << "\n\n";
        // std::cout << num << " x singleSimulationStep           " <<  timeFuncInvocation(num, [&]{singleSimulationStep();}) << " s\n";
        // std::cout << "multipleSimulationSteps(" << num << ")         " <<  timeFuncInvocation(1, [&]{multipleSimulationSteps(num);}) << " s\n";
    }



    void System::cellStep(const Cell& cell)
    {
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
                last_energy_value = totalEnergy(particle);

                if(particle.trySetPosition(particle.position+translation))
                {
                    // energy_after = cell.potentialEnergy(particle, interaction_cutoff);
                    energy_after = totalEnergy(particle);
                                  
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
                energy_after = totalEnergy(particle);
                // energy_after = particle.potentialEnergy(box, interaction_cutoff) 
                //              + interaction::surface_potential(particle, box, interaction_surface_width, interaction_cutoff) * (interaction_surface ? 1 : 0);
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
    }



    Cell& System::cellOfParticle(const Particle& p)
    {
        if(auto cell_it = std::find_if(std::begin(cells), std::end(cells), [&](Cell& c){ return c.contains(p); }); cell_it != std::end(cells))
        {
            return *cell_it;
        }
        else
        {
            throw std::runtime_error("System::cellOfParticle Particle is in no cell");
        }
        
    }



    ParticleRefContainer System::randomParticles(std::size_t num)
    {
        if(num > particles.size())
        {
            throw std::runtime_error("cant find more random particles than the total number of particles");
        }

        ParticleRefContainer sample;
        std::sample(std::begin(particles), std::end(particles), std::back_inserter(sample), num, RandomEngine.pseudo_engine);
        return sample;
    }



    void System::exchangeParticles(Particle& a, Particle& b)
    {
        // std::cout <<"before " << box.distance(a.position, b.position) << "\n";
        b.position = std::exchange(a.position, b.position);
        b.orientation = std::exchange(a.orientation, b.orientation);
        b.sigma = std::exchange(a.sigma, b.sigma);
        b.epsilon = std::exchange(a.epsilon, b.epsilon);
        b.kappa = std::exchange(a.kappa, b.kappa);
        b.gamma = std::exchange(a.gamma, b.gamma);

        b.position_bound_radius_squared = std::exchange(a.position_bound_radius_squared, b.position_bound_radius_squared);
        b.orientation_bound_radiant = std::exchange(a.orientation_bound_radiant, b.orientation_bound_radiant);

        b.surface_affinity_translation = std::exchange(a.surface_affinity_translation, b.surface_affinity_translation);
        b.surface_affinity_rotation = std::exchange(a.surface_affinity_rotation, b.surface_affinity_rotation);

        b.name = std::exchange(a.name, b.name);
        // std::cout <<"after  "  << box.distance(a.position, b.position) << "\n\n";
    }



    void System::globalExchange()
    {
        if(std::isnan(global_exchange_ratio))
        {
            throw std::runtime_error("System::global_exchange_ratio is NaN");
        }
        else if(global_exchange_ratio > 1e-5 && std::isnan(global_exchange_epot_theshold))
        {
            throw std::runtime_error("System::global_exchange_epot_theshold is NaN");
        }

        auto is_valid = [&](const Particle& p) -> bool { 
            return totalEnergy(p)  <  global_exchange_epot_theshold  &&  p.position_bound_radius_squared > std::numeric_limits<REAL>::max()/2  &&  p.orientation_bound_radiant > std::numeric_limits<REAL>::max()/2;
        };

        std::size_t num_particle_exchanges = particles.size()*global_exchange_ratio;
        // if(global_exchange_ratio > 1e-5)
        //     std::cout << "\n\n\n" <<particles.size() << "  " << global_exchange_ratio << "  " << num_particle_exchanges << "\n";

        while(num_particle_exchanges--)
        {
            ParticleRefContainer candidates = randomParticles(2);
            // std::cout << "\nnum_particle_exchanges " << num_particle_exchanges << "\n";
            
            {
                std::size_t not_found = 0; 
                while(!is_valid(candidates.at(0)))
                {
                    candidates[0] = randomParticles(1).at(0);
                    if(++not_found > particles.size())
                    {
                        return;
                    }
                }
                // std::cout << "1 found after " << not_found << "  " << candidates[0].get().repr() << "\n";
            }
            
            {
                std::size_t not_found = 0;
                while((!is_valid(candidates.at(1))) or (candidates[0].get() == candidates[1].get()))
                {
                    candidates[1] = randomParticles(1).at(0);
                    if(++not_found > particles.size())
                    {
                        return;
                    }
                }
                // std::cout << "2 found after " << not_found  << "  " << candidates[1].get().repr() << "\n";
            }
            
            {
                const auto epot_before = totalEnergy(candidates[0]) + totalEnergy(candidates[1]);
                // std::cout << "pre exchange" << "\n" << candidates[0].get().repr() << "\n" << candidates[1].get().repr() << "\n";
                exchangeParticles(candidates[0], candidates[1]);
                // std::cout << "post exchange" << "\n" << candidates[0].get().repr() << "\n" << candidates[1].get().repr() << "\n";
                const auto epot_after = totalEnergy(candidates[0]) + totalEnergy(candidates[1]);

                const bool accepted = ((epot_after - epot_before) < 0) ? true : acceptByMetropolis(epot_after - epot_before, temperature);
                if(!accepted)
                {
                    // std::cout << "exchange declined. delta E " << epot_after - epot_before << " with distance " << box.distance(candidates[0].get().position, candidates[1].get().position) << "\n";
                    exchangeParticles(candidates[0], candidates[1]);
                }
                else
                {
                    // std::cout << "exchange accepted. delta E " << epot_after - epot_before << " with distance " << box.distance(candidates[0].get().position, candidates[1].get().position) << "\n";
                }
            }
        }
        // std::cout << "\n\n";
    }



    void System::makeInteractionLookupTable(ParticleContainer unqiues)
    {
        lookup_table.clear();
        lookup_table = interaction::calculateInteractionLookupTable(unqiues);
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
            .def_readwrite("global_exchange_ratio", &System::global_exchange_ratio)
            .def_readwrite("global_exchange_epot_theshold", &System::global_exchange_epot_theshold)
            .def_readwrite("box", &System::box)
            .def_readwrite("particles", &System::particles)
            .def_readwrite("interaction_cutoff", &System::interaction_cutoff)
            .def_readwrite("interaction_surface", &System::interaction_surface)
            .def_readwrite("interaction_surface_width", &System::interaction_surface_width)
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
            .def("exchangeParticles", &System::exchangeParticles)
            .def("globalExchange", &System::globalExchange)
            .def("randomParticles", &System::randomParticles)
            .def("makeLookupTableFrom", &System::makeInteractionLookupTable)
            .def("makeNeighborLists", &System::makeNeighborLists)
            .def("totalEnergy", static_cast<REAL (System::*)(const Particle&) const>(&System::totalEnergy))
            .def("totalEnergy", static_cast<REAL (System::*)(void)>(&System::totalEnergy))
            .def("potentialEnergy", &System::potentialEnergy)
            .def("potentialEnergyConcurrent", &System::potentialEnergyConcurrent)
            .def("benchmark", &System::benchmark)
            .def("particleEnergies", &System::particleEnergies)
            .def("particleChiValues", &System::particleChiValues)
            .def("particleSurfacePotentialValues", &System::particleSurfacePotentialValues)
            .def("particleExternalPotentialValues", &System::particleExternalPotentialValues)
        ;
    }
}