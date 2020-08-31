#pragma once

#include "definitions.hpp"
#include "type_name.hpp"
#include <Eigen/Core>

#include <pybind11/pybind11.h>
// #include <pybind11/operators.h>

// #pragma GCC diagnostic ignored "-Wdeprecated-declarations" // hurts, but is necessary
// #include <pybind11/eigen.h>
// #pragma GCC diagnostic pop

namespace py = pybind11;


namespace _pyves
{
    struct Parameters
    {
        REAL sys_temperature;
    };



    inline void bind_parameters(py::module& m)
    {
        py::class_<Parameters>(m, "Parameters", py::dynamic_attr())
            .def(py::init<>())
            .def_readwrite("sys_temperature", &Parameters::sys_temperature)
            ;
    }
}