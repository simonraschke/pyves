#pragma once

#include "definitions.hpp"
#include "box.hpp"
#include "observer_ptr.hpp"
#include <atomic>
#include <memory>
#include <Eigen/Geometry>

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>



namespace _pyves 
{ 
    enum class CellState{ UNDEFINED, IDLE, BLOCKED, FINISHED };

    static const CARTESIAN CellBoundOffset
    {
        static_cast<REAL>(1e-3),
        static_cast<REAL>(1e-3),
        static_cast<REAL>(1e-3)
    };

    struct Cell; 
    typedef std::vector<std::reference_wrapper<Cell>> CellRefContainer;
}



PYBIND11_MAKE_OPAQUE(_pyves::CellRefContainer)



struct _pyves::Cell
{
    std::atomic<CellState> state {CellState::UNDEFINED};
    Eigen::AlignedBox<REAL,3> bounding_box;

    CellRefContainer proximity;



    Cell(CARTESIAN_CREF min, CARTESIAN_CREF max)
        // :  bounding_box(std::make_unique<decltype(bounding_box)::element_type>(min-CellBoundOffset, max+CellBoundOffset))
        :  bounding_box(min-CellBoundOffset, max+CellBoundOffset)
    {
    }



    Cell(const Cell &other)  
    { 
        state = other.state.load();
        bounding_box = other.bounding_box;
        proximity = other.proximity;
    }



    Cell& operator=(const Cell &other)  
    { 
        state = other.state.load();
        bounding_box = other.bounding_box;
        proximity = other.proximity;
        return *this; 
    }



    inline bool operator==(const Cell& other) const 
    { 
        return std::addressof(*this) == std::addressof(other); 
    }



    inline bool isNeighbourOf(const Cell& other, const Box<PBC::ON>& b) const
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



    // inline void addToProximity(Cell& other)
    // {
    //     proximity.push_back(std::ref(other));
    // }



    inline bool assertIntegrity() const
    {
        return proximity.size() == 26;
    }



    inline std::string repr() const
    {
        return std::string("<Cell from (") +
            std::to_string(bounding_box.min()(0)) + "|" + std::to_string(bounding_box.min()(1)) + "|" + std::to_string(bounding_box.min()(2)) + 
            ") to (" + 
            std::to_string(bounding_box.max()(0)) + "|" + std::to_string(bounding_box.max()(1)) + "|" + std::to_string(bounding_box.max()(2)) + ")>" ;
    }
};



namespace _pyves
{
    inline void bind_cell(py::module& m)
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
            // .def_property_readonly("bounding_box", [](const Cell& c){ return c.bounding_box.get(); }, py::return_value_policy::reference_internal)
            // .def_readonly("proximity", &Cell::proximity)            
            // .def_property("proximity", &Cell::proximity, [](Cell& i, Cell& j ) { i.addToProximity(j); }, py::return_value_policy::reference_internal)
            .def_readwrite("proximity", &Cell::proximity, py::return_value_policy::reference_internal)
            .def_property("min", [](const Cell& c){ return c.bounding_box.min(); }, [](Cell& c, CARTESIAN_CREF v) { c.bounding_box.min() = v - CellBoundOffset; }, py::return_value_policy::reference_internal)
            .def_property("max", [](const Cell& c){ return c.bounding_box.max(); }, [](Cell& c, CARTESIAN_CREF v) { c.bounding_box.max() = v - CellBoundOffset; }, py::return_value_policy::reference_internal)
            .def("isNeighbourOf", &Cell::isNeighbourOf)
            .def("assertIntegrity", &Cell::assertIntegrity)
            .def("__repr__", &Cell::repr)
        ;
    }
} // namespace _pyves