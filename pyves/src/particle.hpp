#pragma once

#include "definitions.hpp"
#include "type_name.hpp"
#include <Eigen/Core>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#pragma GCC diagnostic ignored "-Wdeprecated-declarations" // hurts, but is necessary
#include <pybind11/eigen.h>
#pragma GCC diagnostic pop

namespace py = pybind11;


namespace _pyves
{
    struct Particle
    {
        CARTESIAN position;
        CARTESIAN orientation;
        std::string name = "UNDEF";
        

        inline auto x() -> decltype(position.coeffRef(0)) { return position.coeffRef(0); }
        inline auto y() -> decltype(position.coeffRef(1)) { return position.coeffRef(1); }
        inline auto z() -> decltype(position.coeffRef(2)) { return position.coeffRef(2); }     

        inline auto ux() -> decltype(orientation.coeffRef(0)) { return orientation.coeffRef(0); }
        inline auto uy() -> decltype(orientation.coeffRef(1)) { return orientation.coeffRef(1); }
        inline auto uz() -> decltype(orientation.coeffRef(2)) { return orientation.coeffRef(2); }

        inline bool operator==(const Particle& other) const { return std::addressof(*this) == std::addressof(other); };
        inline bool operator!=(const Particle& other) const { return std::addressof(*this) != std::addressof(other); };
        
        inline std::string repr()
        {
            return 
                std::string("<Particle (REAL=") + 
                type_name<REAL>() + 
                ") at (" + 
                std::to_string(x()) + "|" + std::to_string(y()) + "|" + std::to_string(z()) + 
                ") in (" + 
                std::to_string(ux()) + "|" + std::to_string(uy()) + "|" + std::to_string(uz()) + 
                ") direction>";
        }  
    };



    inline void bind_particle(py::module& m) 
    {
        py::class_<Particle>(m, "Particle", py::dynamic_attr())
            .def(py::init<>())
            .def(py::init<Particle>())
            .def(py::init<CARTESIAN_CREF, CARTESIAN_CREF>())
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def_readwrite("position", &Particle::position)
            .def_readwrite("orientation", &Particle::orientation)
            .def_readwrite("name", &Particle::name)
            .def("x", &Particle::x, py::return_value_policy::reference_internal)
            .def("y", &Particle::y, py::return_value_policy::reference_internal)
            .def("z", &Particle::z, py::return_value_policy::reference_internal)
            .def("__repr__", &Particle::repr);
            ;
    }
}