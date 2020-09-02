#pragma once

#include "box.hpp"
#include "parameters.hpp"



namespace _pyves
{
    struct System
    {
        Box<PBC::ON> box;
        Parameters prms;
        // particles
        // cells
        // 
    };



    inline void bind_system(py::module& m)
    {
        py::class_<System>(m, "System", py::dynamic_attr())
            .def(py::init<>())
            .def_readwrite("prms", &System::prms)
            ;
    }
}