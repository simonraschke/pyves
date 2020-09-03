#pragma once

#include "definitions.hpp"
// #include <deque>
// #include <numeric>
#include <atomic>
#include <vector>
#include <Eigen/Geometry>
// #include <tbb/spin_rw_mutex.h>
// #include <tbb/concurrent_vector.h>

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



PYBIND11_MAKE_OPAQUE(std::vector<std::reference_wrapper<_pyves::Cell>>)



struct _pyves::Cell
{
    std::atomic<CellState> state {CellState::UNDEFINED};
    Eigen::AlignedBox<REAL,3> bounding_box;

    CellRefContainer proximity;

    Cell(CARTESIAN_CREF from, CARTESIAN_CREF to)
        :  bounding_box(from-CellBoundOffset, to+CellBoundOffset)
    {
    }

    std::string repr() const
    {
        return std::string("<Cell>");
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
        
        py::class_<Cell>(m, "Cell")
            .def(py::init<CARTESIAN_CREF,CARTESIAN_CREF>(), py::arg("from"), py::arg("to"))
            .def_readonly("bounding_box", &Cell::bounding_box)
            .def_readonly("proximity", &Cell::proximity)
            ;
    }
} // namespace _pyves