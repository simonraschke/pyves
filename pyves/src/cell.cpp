#include "cell.hpp"

namespace _pyves
{    
    Cell::Cell(CARTESIAN_CREF min, CARTESIAN_CREF max)
        // :  bounding_box(std::make_unique<decltype(bounding_box)::element_type>(min-CellBoundOffset, max+CellBoundOffset))
        :  bounding_box(min-CellBoundOffset, max+CellBoundOffset)
    {
    }



    Cell::Cell(const Cell &other)  
    { 
        state = other.state.load();
        bounding_box = other.bounding_box;
        proximity = other.proximity;
    }



    Cell& Cell::operator=(const Cell &other)  
    { 
        state = other.state.load();
        bounding_box = other.bounding_box;
        proximity = other.proximity;
        return *this; 
    }



    bool Cell::operator==(const Cell& other) const 
    { 
        return std::addressof(*this) == std::addressof(other); 
    }



    bool Cell::isNeighbourOf(const Cell& other, const Box<PBC::ON>& b) const
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
            const CARTESIAN connection_vector = b.distanceVector(bounding_box.center(), other.bounding_box.center()).cwiseAbs();
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



    bool Cell::contains(const Particle& p) const
    {
        return bounding_box.contains(p.position);
    }



    bool Cell::assertIntegrity() const
    {
        return all(
            proximity.size() == 26,
            region.size() == 27
        );
    }



    void Cell::shuffle()
    {
        std::shuffle(std::begin(particles), std::end(particles), std::mt19937_64(std::random_device{}()));
        std::shuffle(std::begin(proximity), std::end(proximity), std::mt19937_64(std::random_device{}()));
        std::shuffle(std::begin(region), std::end(region), std::mt19937_64(std::random_device{}()));
    }



    REAL Cell::potentialEnergy(const Particle& particle, const Box<PBC::ON>& box, REAL cutoff) const
    {
        return std::accumulate(std::begin(region), std::end(region), REAL(0), [&](REAL __val, const Cell& cell)
        {
            return __val + std::accumulate(std::begin(cell.particles), std::end(cell.particles), REAL(0), [&](REAL _val, const Particle& compare)
            {
                return particle == compare ? _val : _val + interaction(particle, compare, box, cutoff);
            });
        });
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
            .def(py::init<CARTESIAN_CREF,CARTESIAN_CREF>(), py::arg("min"), py::arg("max"))
            .def_readonly("bounding_box", &Cell::bounding_box)
            .def_readwrite("proximity", &Cell::proximity, py::return_value_policy::reference_internal)
            .def_readwrite("region", &Cell::region, py::return_value_policy::reference_internal)
            .def_readwrite("particles", &Cell::particles, py::return_value_policy::reference_internal)
            .def_property("min", [](const Cell& c){ return c.bounding_box.min(); }, 
                                 [](Cell& c, CARTESIAN_CREF v) { c.bounding_box.min() = v - CellBoundOffset; }, py::return_value_policy::reference_internal)
            .def_property("max", [](const Cell& c){ return c.bounding_box.max(); }, 
                                 [](Cell& c, CARTESIAN_CREF v) { c.bounding_box.max() = v - CellBoundOffset; }, py::return_value_policy::reference_internal)
            .def("contains", &Cell::contains)
            .def("isNeighbourOf", &Cell::isNeighbourOf)
            .def("assertIntegrity", &Cell::assertIntegrity)
            .def("shuffle", &Cell::shuffle)
            .def("__repr__", &Cell::repr)
        ;
    }
}