#include "system.hpp"


namespace _pyves
{
    void System::setThreads(std::size_t num)
    {
        threads = num;
        task_arena.initialize(threads);
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



    bool System::assertIntegrity() const
    {
        return all(
            std::all_of(std::begin(particles), std::end(particles), [](const auto& p) { return p.assertIntegrity(); }),
            std::all_of(std::begin(cells), std::end(cells), [](const auto& c) { return c.assertIntegrity(); }),
            std::isfinite(cores),
            std::isfinite(threads),
            std::isfinite(temperature),
            std::isfinite(time_max)
        );
    }



    void System::singleSimulationStep()
    {    
        cellBasedApplyFunctor([&](const auto& c){ cellStep(c);});
    }



    void System::cellStep(const Cell& cell)
    {
        // TODO:FIXME:TODO:FIXME: CAUTION: this function is exactly what it should look like.
        // do not change order or modify any function call

        REAL last_energy_value;
        REAL energy_after;
        std::uniform_real_distribution<REAL> dist_coords(-translation_alignment(),translation_alignment());
        std::uniform_real_distribution<REAL> dist_orientation(-rotation_alignment(),rotation_alignment());
        CARTESIAN translation;
        CARTESIAN orientation_before;

        std::vector<ParticleRefContainer::const_iterator> iterators(cell.particles.size());
        std::iota(std::begin(iterators), std::end(iterators), std::begin(cell.particles));
        std::shuffle(std::begin(iterators), std::end(iterators), pseudo_engine);

        for(const auto& iterator : iterators)
        {
            Particle& particle = *iterator;
            
            // coordinates move
            {
                translation = CARTESIAN
                (
                    dist_coords(pseudo_engine),
                    dist_coords(pseudo_engine),
                    dist_coords(pseudo_engine)
                );

                last_energy_value = cell.potentialEnergy(particle, box, 3);
                // if(particle->try_setCoordinates(particle->getCoordinates()+translation))
                if(particle.position += translation; true)
                {
                    energy_after = cell.potentialEnergy(particle, box, 3);

                    // rejection
                    if(!acceptByMetropolis(energy_after - last_energy_value, temperature))
                    {
                        translation_alignment.rejected();
                        particle.position -= translation;
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

            CARTESIAN random_vector = CARTESIAN::Random();

            // to not oversample with points outside of the unity sphere
            while(random_vector.squaredNorm() > 1)
            {
                random_vector = CARTESIAN::Random();
            }

            // if( 
            //     orientation_before = particle.orientation; 
            //     particle->try_setOrientation(Eigen::AngleAxis<REAL>(dist_orientation(pseudo_engine), random_vector) * particle->getOrientation())
            // )
            // {
            //     energy_after = cell.potentialOfSingleParticle(*particle);

            //     // rejection
            //     if(!acceptance.isValid(energy_after - last_energy_value))
            //     {
            //         sw_orientation.rejected();
            //         particle->getOrientation() = orientation_before;
            //     }
            //     // acctance
            //     else
            //     {
            //         sw_orientation.accepted();
            //         last_energy_value = energy_after;
            //     }
            // }
            // else
            // {
            //     sw_orientation.rejected();
            // }
        }

        // TODO:FIXME:TODO:FIXME: CAUTION: this function is exactly what it should look like.
        // do not change order or modify any function call
    }



    


    void System::shuffle()
    {
        tbb::task_group tg;
        task_arena.execute([&]()
        {
            tg.run_and_wait([&]{
                tbb::parallel_for_each(std::begin(cells), std::end(cells), [&] (Cell& cell) { cell.shuffle(); });
            });
        });
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
            .def_readwrite("particles", &System::particles, py::return_value_policy::reference_internal)
            .def_readwrite("cells", &System::cells, py::return_value_policy::reference_internal)
            .def("particleIsFree", static_cast<bool (System::*)(const Particle&) const>(&System::particleIsFree))
            .def("particleIsFree", static_cast<bool (System::*)(const Particle&, REAL) const>(&System::particleIsFree))
            .def("shuffle", &System::shuffle)
            .def_readwrite("translation", &System::translation_alignment, py::return_value_policy::reference_internal)
            .def_readwrite("rotation", &System::rotation_alignment, py::return_value_policy::reference_internal)
            .def("assertIntegrity", &System::assertIntegrity)
            .def("singleSimulationStep", &System::singleSimulationStep)
        ;
    }
}