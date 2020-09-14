#include "system.hpp"


namespace _pyves
{
    void System::setThreads(std::size_t num)
    {
        // threads = std::max(num, static_cast<std::size_t>(2));
        threads = num;
	    // tbb::global_control control(tbb::global_control::max_allowed_parallelism, threads);

        // if(threads == 1)
        // {
        //     task_arena.initialize(threads, 1);
        // }
        // else
        // {
        //     task_arena.initialize(threads);
        // }
        executor.reset(new tf::Executor(threads));
    }



    bool System::particleIsFree(const Particle& subject) const
    {
        return std::all_of(std::begin(particles), std::end(particles), [&](const auto& p)
        {
            if(std::isnan(subject.sigma) || std::isnan(p.sigma))
            {
                throw std::logic_error("System::particleIsFree found particle where sigma is NaN");
            }
            return subject == p ? true : box.distanceParticle(subject, p) > 1.1224f*(subject.sigma+p.sigma)/2;
        });
    }



    bool System::particleIsFree(const Particle& subject, REAL cutoff) const
    {
        const REAL cutoff_squared = cutoff*cutoff;
        return std::all_of(std::begin(particles), std::end(particles), [&](const auto& p)
        { 
            return subject == p ? true : box.squaredDistanceParticle(subject, p) > cutoff_squared;
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
            std::isfinite(cores),
            std::isfinite(threads),
            std::isfinite(temperature)
        );
    }



    void System::singleSimulationStep()
    {
        prepareSimulationStep();
        applyToCells([&](const Cell& c){ cellStep(c);});
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
        tf::Taskflow taskflow;
        taskflow.for_each_static(std::begin(cells), std::end(cells), [&] (Cell& cell) 
        {
            cell.shuffle(); 
        });
        executor->run(taskflow).get();
    }



    void System::prepareSimulationStep()
    {
        // tf::Executor executor(threads);
        tf::Taskflow taskflow;
        taskflow.for_each_static(std::begin(cells), std::end(cells), [&] (Cell& cell) 
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
            cell.state = (cell.particles.size() > 0) ? CellState::IDLE : CellState::FINISHED;
            cell.shuffle(); 
        });

        executor->run(taskflow).get();
    }



    void System::cellStep(const Cell& cell)
    {
        // TODO:FIXME:TODO:FIXME: CAUTION: this function is exactly what it should look like.
        // do not change order or modify any function call

        // std::cout << cell.repr() << "\n";
        REAL last_energy_value;
        REAL energy_after;
        // std::uniform_real_distribution<REAL> dist_coords(-translation_alignment(),translation_alignment());
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
                    translation = CARTESIAN::Random() * translation_alignment();
                }
                while(translation.squaredNorm() > translation_alignment()*translation_alignment());

                last_energy_value = cell.potentialEnergy(particle, 3);
                
                if(particle.trySetPosition(particle.position+translation))
                {
                    energy_after = cell.potentialEnergy(particle, 3);

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
                energy_after = cell.potentialEnergy(particle, 3);

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



        py::class_<System>(m, "System", py::dynamic_attr())
            .def(py::init<>())
            .def_readwrite("cores", &System::cores, py::return_value_policy::reference_internal)
            .def_property("threads",[](const System& s){ return s.threads; }, &System::setThreads)
            .def_readwrite("temperature", &System::temperature, py::return_value_policy::reference_internal)
            .def_readwrite("box", &System::box, py::return_value_policy::reference_internal)
            .def_readwrite("particles", &System::particles)
            .def_readwrite("cells", &System::cells, py::return_value_policy::reference_internal)
            .def("particleIsFree", static_cast<bool (System::*)(const Particle&) const>(&System::particleIsFree))
            .def("particleIsFree", static_cast<bool (System::*)(const Particle&, REAL) const>(&System::particleIsFree))
            .def("shuffle", &System::shuffle)
            .def("numParticlesInCells", &System::numParticlesInCells)
            .def_readwrite("translation", &System::translation_alignment, py::return_value_policy::reference_internal)
            .def_readwrite("rotation", &System::rotation_alignment, py::return_value_policy::reference_internal)
            .def("assertIntegrity", &System::assertIntegrity)
            .def("prepareSimulationStep", &System::prepareSimulationStep)
            .def("singleSimulationStep", &System::singleSimulationStep)
            .def("multipleSimulationSteps", &System::multipleSimulationSteps)
        ;
    }
}