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
        localExchange();
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
        return p.potentialEnergy(box, interaction_cutoff) + interaction::external_potential(p, box, interaction_surface_width, interaction_cutoff) * (interaction_surface ? 1 : 0);
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



    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> System::distanceMatrix() const
    {
        Eigen::Matrix<REAL,Eigen::Dynamic, Eigen::Dynamic> matrix (particles.size(), particles.size());
        matrix.setZero();
        
        for(std::size_t i = 0; i < particles.size(); ++i)
        {
            for(std::size_t j = 0; j < i; ++j)
            {
                matrix(i,j) = box.distance(particles[i].getPosition(), particles[j].getPosition());
                matrix(j,i) = matrix(i,j);
            }
        }

        return matrix;
    }



    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> System::potentialEnergyMatrix() const
    {
        Eigen::Matrix<REAL,Eigen::Dynamic, Eigen::Dynamic> matrix (particles.size(), particles.size());
        matrix.setZero();

        for(std::size_t i = 0; i < particles.size(); ++i)
        {
            matrix(i,i) = std::numeric_limits<REAL>::quiet_NaN();
            // matrix(i,i) = std::numeric_limits<REAL>::infinity();
        }
        
        for(std::size_t i = 0; i < particles.size(); ++i)
        {
            for(std::size_t j = 0; j < i; ++j)
            {
                matrix(i,j) = interaction::potentialEnergy(particles[i], particles[j], box, interaction_cutoff);
                matrix(j,i) = matrix(i,j);
            }
        }

        return matrix;
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



    Particle& System::randomParticle()
    {
        return random(particles);
    }



    void System::exchangeParticleParameters(Particle& a, Particle& b, bool orientation = false)
    {
        // 2 possible ways:
        // a) either exchange position and neighbour list (too complicated with large overhead combined with linked cell-lists)
        // b) or exchange parameters and name (done here)

        if(a.position_bound_radius_squared < std::numeric_limits<REAL>::max()/2  ||  a.orientation_bound_radiant < std::numeric_limits<REAL>::max()/2)
        {
            throw std::logic_error("cant exchange particle with restricted movement");
        }
        else if(b.position_bound_radius_squared < std::numeric_limits<REAL>::max()/2  ||  b.orientation_bound_radiant < std::numeric_limits<REAL>::max()/2)
        {
            throw std::logic_error("cant exchange particle with restricted movement");
        }
        else if(a == b)
        {
            throw std::logic_error("cant exchange one particle with itself");
        }

        // paramters
        if(orientation)
        {
            b.orientation = std::exchange(a.orientation, b.orientation);
        }

        b.sigma = std::exchange(a.sigma, b.sigma);
        b.epsilon = std::exchange(a.epsilon, b.epsilon);
        b.kappa = std::exchange(a.kappa, b.kappa);
        b.gamma = std::exchange(a.gamma, b.gamma);
        b.surface_affinity_translation = std::exchange(a.surface_affinity_translation, b.surface_affinity_translation);
        b.surface_affinity_rotation = std::exchange(a.surface_affinity_rotation, b.surface_affinity_rotation);
        b.self_affinity = std::exchange(a.self_affinity, b.self_affinity);
        b.other_affinity = std::exchange(a.other_affinity, b.other_affinity);
        
        std::swap(a.initial_position, b.initial_position);
        std::swap(a.initial_orientation, b.initial_orientation);

        // name
        b.name = std::exchange(a.name, b.name);
    }


    
    void System::localExchange()
    {
        static const bool exchange_is_possible = [&]{
            using std::string;
            std::unordered_set<std::string> set;
            for (const Particle& p : particles)
            {
                if(   p.position_bound_radius_squared > std::numeric_limits<REAL>::max()/2  
                   && p.orientation_bound_radiant > std::numeric_limits<REAL>::max()/2
                )
                {
                    set.insert(p.name);
                }
            }
            const string explanation = (set.size() >= 2) ? string("possible") : string("impossible");
            std::cout << "local exchange is " << explanation << " with " << set.size() << " unique exchangeable particles\n";
            return set.size() >= 2;
        }();

        if(!exchange_is_possible)
        {
            return;
        }
        else if( exchange_local_number == 0 )
        {
            return;
        }
        else if((exchange_local_number > 0) && std::isnan(exchange_local_etot_theshold))
        {
            throw std::runtime_error("System::exchange_local_etot_theshold is NaN");
        }


        auto is_valid_compare = [&](const Particle& compare, const Particle& origin) -> bool { 
            return (totalEnergy(compare) < exchange_local_etot_theshold)  
                && (compare.position_bound_radius_squared > std::numeric_limits<REAL>::max()/2 ) 
                && (compare.orientation_bound_radiant > std::numeric_limits<REAL>::max()/2 )
                && (compare.name != origin.name);
        };

        auto is_valid_origin = [&](const Particle& origin) -> bool { 
            return (totalEnergy(origin) < exchange_local_etot_theshold)
                && (origin.position_bound_radius_squared > std::numeric_limits<REAL>::max()/2)
                && (origin.orientation_bound_radiant > std::numeric_limits<REAL>::max()/2) 
                && std::any_of(std::begin(origin.neighbors), std::end(origin.neighbors), [&] (const Particle& n){ return is_valid_compare(n, origin); })
            ;
        };



        for(decltype(exchange_local_number) i = 0; i < exchange_local_number; ++i)
        {
            std::reference_wrapper<Particle> origin = random(particles);

            {
                auto limit = particles.size();
                while(!is_valid_origin(origin))
                {
                    if(! (--limit))
                    {
                        return;
                    }
                    origin = random(particles);
                }
            }

            bool no_neighbor_found = false;
            std::reference_wrapper<Particle> compare = random(origin.get().neighbors);

            {
                auto limit = particles.size();
                while(!is_valid_compare(compare, origin))
                {
                    if(! (--limit))
                    {
                        no_neighbor_found = true;
                        break;
                    }
                    compare = random(origin.get().neighbors);
                }
            }

            if(no_neighbor_found)
            {
                continue;
            }
            else
            {
                const auto epot_before = totalEnergy(origin) + totalEnergy(compare);
                exchangeParticleParameters(origin, compare, exchange_local_orientation);
                const auto epot_after = totalEnergy(origin) + totalEnergy(compare);

                const bool accepted = ((epot_after - epot_before) < 0) ? true : acceptByMetropolis(epot_after - epot_before, temperature);
                if(!accepted)
                {
                    exchangeParticleParameters(origin, compare, exchange_local_orientation);
                }
                else
                {
                }
            }
        }
    }


    
    void System::globalExchange()
    {
        static const bool exchange_is_possible = [&]{
            using std::string;
            std::unordered_set<std::string> set;
            for (const Particle& p : particles)
            {
                if(   p.position_bound_radius_squared > std::numeric_limits<REAL>::max()/2  
                   && p.orientation_bound_radiant > std::numeric_limits<REAL>::max()/2
                )
                {
                    set.insert(p.name);
                }
            }
            const string explanation = (set.size() >= 2) ? string("possible") : string("impossible");
            std::cout << "global exchange is " << explanation << " with " << set.size() << " unique exchangeable particles\n";
            return set.size() >= 2;
        }();

        if(!exchange_is_possible)
        {
            return;
        }
        else if( exchange_global_number == 0 )
        {
            return;
        }
        else if((exchange_global_number > 0) && std::isnan(exchange_global_etot_theshold))
        {
            throw std::runtime_error("System::exchange_global_etot_theshold is NaN");
        }

        auto is_valid_origin = [&](const Particle& origin) -> bool { 
            return (totalEnergy(origin)  <  exchange_global_etot_theshold)  
                && (origin.position_bound_radius_squared > std::numeric_limits<REAL>::max()/2)  
                && (origin.orientation_bound_radiant > std::numeric_limits<REAL>::max()/2);
        };

        auto is_valid_compare = [&](const Particle& compare, const Particle& origin) -> bool { 
            return (totalEnergy(compare)  <  exchange_global_etot_theshold)  
                && (compare.position_bound_radius_squared > std::numeric_limits<REAL>::max()/2)  
                && (compare.orientation_bound_radiant > std::numeric_limits<REAL>::max()/2) 
                && (compare.name != origin.name);
        };



        for(decltype(exchange_global_number) i = 0; i < exchange_global_number; ++i)
        {
            std::reference_wrapper<Particle> origin = randomParticle();
            std::reference_wrapper<Particle> compare = randomParticle();

            {   
                auto limit = particles.size();
                while(!is_valid_origin(origin))
                {
                    if(! (--limit))
                    {
                        return;
                    }
                    origin = randomParticle();
                }
            }

            {
                auto limit = particles.size();
                while(!is_valid_compare(compare, origin) or (compare.get() == origin.get()))
                {
                    if(! (--limit))
                    {
                        return;
                    }
                    compare = randomParticle();
                }
            }

            {
                const auto epot_before = totalEnergy(origin) + totalEnergy(compare);
                exchangeParticleParameters(origin, compare, exchange_global_orientation);
                const auto epot_after = totalEnergy(origin) + totalEnergy(compare);

                const bool accepted = ((epot_after - epot_before) < 0) ? true : acceptByMetropolis(epot_after - epot_before, temperature);
                if(!accepted)
                {
                    exchangeParticleParameters(origin, compare, exchange_global_orientation);
                }
            }
        }
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

        py::class_<System>(m, "System")//, py::dynamic_attr())
            .def(py::init<>())
            .def_property("threads", [](const System& s){ return s.threads; }, &System::setThreads)
            .def_readwrite("temperature", &System::temperature)
            .def_readwrite("exchange_global_number", &System::exchange_global_number)
            .def_readwrite("exchange_global_etot_theshold", &System::exchange_global_etot_theshold)
            .def_readwrite("exchange_global_orientation", &System::exchange_global_orientation)
            .def_readwrite("exchange_local_number", &System::exchange_local_number)
            .def_readwrite("exchange_local_etot_theshold", &System::exchange_local_etot_theshold)
            .def_readwrite("exchange_local_orientation", &System::exchange_local_orientation)
            .def_readwrite("box", &System::box)
            .def_readwrite("particles", &System::particles)
            .def_readwrite("interaction_cutoff", &System::interaction_cutoff)
            .def_readwrite("interaction_surface", &System::interaction_surface)
            .def_readwrite("interaction_surface_width", &System::interaction_surface_width)
            .def_readwrite("cells", &System::cells)
            .def_readwrite("cell_update_interval", &System::cell_update_interval)
            .def_readwrite("neighbor_cutoff", &System::neighbor_cutoff)
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
            .def("exchangeParticleParameters", &System::exchangeParticleParameters, py::arg("a"), py::arg("b"), py::arg("orientation") = bool(false))
            .def("localExchange", &System::localExchange)
            .def("globalExchange", &System::globalExchange)
            .def("randomParticle", &System::randomParticle, py::return_value_policy::reference)
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
            .def("distanceMatrix", &System::distanceMatrix)
            .def("potentialEnergyMatrix", &System::potentialEnergyMatrix)
        ;
    }
}