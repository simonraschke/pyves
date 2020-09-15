#include "cell.hpp"

namespace _pyves
{    
    Cell::Cell(CARTESIAN_CREF min, CARTESIAN_CREF max, const Box<PBC::ON>& _box)
        // :  bounding_box(std::make_unique<decltype(bounding_box)::element_type>(min-CellBoundOffset, max+CellBoundOffset))
        : bounding_box(min-CellBoundOffset, max+CellBoundOffset)
        , box(_box)
    {
    }



    Cell::Cell(const Cell& other)
    { 
        state = other.state.load();
        bounding_box = other.bounding_box;
        proximity = other.proximity;
        box = other.box;
    }



    Cell& Cell::operator=(const Cell &other)  
    { 
        state = other.state.load();
        bounding_box = other.bounding_box;
        proximity = other.proximity;
        box = other.box;
        return *this; 
    }



    bool Cell::operator==(const Cell& other) const 
    { 
        return this == &other; 
    }



    void Cell::removeParticle(const Particle& to_remove)
    {
        std::lock_guard<std::shared_mutex> lock(particles_access_mutex);
        // tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, true);
        particles.erase( std::remove_if(std::begin(particles), std::end(particles), [&](const Particle& to_compare)
        { 
            return to_remove == to_compare;
        ;}), std::end(particles));

    }
    


    bool Cell::try_add(Particle& particle)
    {
        if(!contains(particle) && insideCellBounds(box.scaleToBox(CARTESIAN(particle.position))))
        {   
            // tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, true);
            // particles.push_back(std::ref(particle));
            // lock.release();
            {
                std::lock_guard<std::shared_mutex> lock(particles_access_mutex);
                particles.push_back(std::ref(particle));
            }
            // std::cout << "   "  << __func__ << "  Cell contains particle " << box.scaleToBox(CARTESIAN(particle.position)).format(VECTORFORMAT) << " after insertion " << std::boolalpha << contains(particle) << "\n";
            assert(contains(particle));
            return true;
        }
        else
        {
            // std::cout << "   "   << __func__ << "  Cell contains particle " << box.scaleToBox(CARTESIAN(particle.position)).format(VECTORFORMAT) << " after NO insertion " << std::boolalpha << contains(particle) << "\n";
            assert(!contains(particle));
            return false;
        }
    }



    bool Cell::isNeighbourOf(const Cell& other) const
    {
        if(*this == other)
        {
            return false;
        }
        else if(bounding_box.intersects(other.bounding_box))
        {
            return true;
        }
        else
        {
            // std::cout << box.getLengthX() << " "<<  box.getLengthY() << " "<<  box.getLengthZ() << "\n";
            const CARTESIAN connection_vector = box.distanceVector(bounding_box.center(), other.bounding_box.center()).cwiseAbs();
            if(     (connection_vector(0) < bounding_box.sizes()(0) + 0.001) 
                &&  (connection_vector(1) < bounding_box.sizes()(1) + 0.001) 
                &&  (connection_vector(2) < bounding_box.sizes()(2) + 0.001))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }



    bool Cell::insideCellBounds(const Particle& p) const
    {
        return bounding_box.contains(box.scaleToBox(p.position));
    }



    bool Cell::insideCellBounds(CARTESIAN_CREF c) const
    {
        return bounding_box.contains(c);
    }



    bool Cell::contains(const Particle& p)
    {
        // tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, false);
        std::shared_lock<std::shared_mutex> lock(particles_access_mutex);
        return std::any_of(std::begin(particles), std::end(particles), [&](const Particle& other ){ return other == p; });
    }



    bool Cell::assertIntegrity()
    {
        // std::shared_lock<std::shared_mutex> lock(particles_access_mutex);
        return all(
            proximity.size() == 26,
            region.size() == 27,
            std::all_of(std::begin(particles), std::end(particles), [&](const Particle& p){ return contains(p); }),
            std::all_of(std::begin(particles), std::end(particles), [&](const Particle& p){ return insideCellBounds(p); })
        );
    }



    void Cell::shuffle()
    {
        // tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, true);
        std::lock_guard<std::shared_mutex> lock(particles_access_mutex);
        std::shuffle(std::begin(particles), std::end(particles), RandomEngine.pseudo_engine);
        std::shuffle(std::begin(proximity), std::end(proximity), RandomEngine.pseudo_engine);
        std::shuffle(std::begin(region), std::end(region), RandomEngine.pseudo_engine);
    }



    auto Cell::particlesOutOfBounds() -> std::deque<decltype(particles)::value_type>
    {
        // tbb::spin_rw_mutex::scoped_lock lock(particles_access_mutex, false);
        std::shared_lock<std::shared_mutex> lock(particles_access_mutex);
        std::deque<decltype(particles)::value_type> leavers;
        for(Particle& particle : particles)
        {
            if(!insideCellBounds(box.scaleToBox(particle.position)))
            {
                leavers.push_back(std::ref(particle));
            }
        }
        return leavers;
    }



    REAL Cell::potentialEnergy(const Particle& particle, REAL cutoff) const
    {
        return std::accumulate(std::begin(region), std::end(region), REAL(0), [&](REAL __val, const Cell& cell)
        {
            return __val + std::accumulate(std::begin(cell.particles), std::end(cell.particles), REAL(0), [&](REAL _val, const Particle& compare)
            {
                return particle == compare ? _val : _val + interaction(particle, compare, box, cutoff);
            });
        });
    }



    REAL Cell::potentialEnergyVectorized(const Particle& particle, REAL cutoff) const
    {
        // const std::size_t num_particles = std::accumulate(std::begin(region), std::end(region), std::size_t(0), [](auto i, const Cell& cell){ return i + cell.particles.size(); });
        
        // decltype(Cell::particles) region_particles;
        // region_particles.reserve(num_particles);

        // for(Cell& c : region)
        // {
        //     std::transform(c.particles.begin(), c.particles.end(), std::back_inserter(region_particles), [](auto& x){ return std::ref<Particle>(x); });
        // }
        
        // std::cout << "\n";
        // std::cout << "\n";
        // std::cout << "\n";
        // Eigen::MatrixXf positions(num_particles, 3);
        // for(std::size_t i = 0; i < region_particles.size(); ++i)
        // {
        //     positions(i,0) = region_particles.at(i).get().getx();
        //     positions(i,1) = region_particles.at(i).get().gety();
        //     positions(i,2) = region_particles.at(i).get().getz();
        // }
        // std::cout << positions.format(PYTHONFORMAT) << "\n";

        // auto distance_vectors = positions.rowwise() - particle.position.transpose();
        // std::cout << distance_vectors.format(PYTHONFORMAT) << "\n";

        return 0;
    }



    std::string Cell::repr() const
    {
        std::stringstream ss;
        ss  << "<Cell from " << bounding_box.min().format(VECTORFORMAT) << " to " << bounding_box.max().format(VECTORFORMAT) << ">";
        return ss.str();
    }





    void bind_cell(py::module& m)
    {
        py::enum_<CellState>(m, "CellState")
            .value("UNDEFINED", CellState::UNDEFINED)
            .value("IDLE", CellState::IDLE)
            .value("BLOCKED", CellState::BLOCKED)
            .value("FINISHED", CellState::FINISHED)
            .export_values()
        ;



        py::bind_vector<CellRefContainer>( m, "CellRefContainer" )
            .def(py::init<>())
            .def("__len__", [](const CellRefContainer& v) { return v.size(); })
            .def("__iter__", [](CellRefContainer& v) 
            {
                return py::make_iterator(std::begin(v), std::end(v));
            }, py::keep_alive<0, 1>())
            .def("__repr__", [](const CellRefContainer& v) {
                return "CellRefContainer\n[\n"
                    + std::accumulate(std::begin(v), std::end(v), std::string(""), [](std::string s, const Cell& c) 
                    { 
                        return s+"\t"+c.repr()+",\n";
                    })
                    + "]";
            })
        ;



        m.def("CellBoundOffset", [](){ return CellBoundOffset; });



        using BoxClass = Eigen::AlignedBox<REAL,3>;
        py::class_<BoxClass>(m, "AlignedBox")
            .def(py::init<>())
            .def(py::init<CARTESIAN_CREF,CARTESIAN_CREF>())
            .def_property_readonly("center", &BoxClass::center, py::return_value_policy::reference_internal)
            .def_property("min", [](const BoxClass& b){ return b.min(); }, [](BoxClass& b, CARTESIAN_CREF c) { b.min() = c; }, py::return_value_policy::reference_internal)
            .def_property("max", [](const BoxClass& b){ return b.max(); }, [](BoxClass& b, CARTESIAN_CREF c) { b.max() = c; }, py::return_value_policy::reference_internal)
        ;



        py::class_<Cell>(m, "Cell", py::dynamic_attr())
            .def(py::init<CARTESIAN_CREF, CARTESIAN_CREF, Box<PBC::ON>>(), py::arg("min"), py::arg("max"), py::arg("box"))
            .def_readonly("bounding_box", &Cell::bounding_box)
            .def_readwrite("proximity", &Cell::proximity, py::return_value_policy::reference_internal)
            .def_readwrite("region", &Cell::region, py::return_value_policy::reference_internal)
            .def_readwrite("particles", &Cell::particles, py::return_value_policy::reference_internal)
            .def_property("min", [](const Cell& c){ return c.bounding_box.min(); }, 
                                 [](Cell& c, CARTESIAN_CREF v) { c.bounding_box.min() = v - CellBoundOffset; }, py::return_value_policy::reference_internal)
            .def_property("max", [](const Cell& c){ return c.bounding_box.max(); }, 
                                 [](Cell& c, CARTESIAN_CREF v) { c.bounding_box.max() = v - CellBoundOffset; }, py::return_value_policy::reference_internal)
            .def("particlesOutOfBounds", &Cell::particlesOutOfBounds)
            .def("insideCellBounds", static_cast<bool (Cell::*)(const Particle&) const>(&Cell::insideCellBounds))
            .def("insideCellBounds", static_cast<bool (Cell::*)(CARTESIAN_CREF) const>(&Cell::insideCellBounds))
            .def("contains", &Cell::contains)
            .def("isNeighbourOf", &Cell::isNeighbourOf)
            .def("assertIntegrity", &Cell::assertIntegrity)
            .def("shuffle", &Cell::shuffle)
            .def("__repr__", &Cell::repr)
        ;
    }
}