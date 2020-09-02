#pragma once

#include "definitions.hpp"
#include "type_name.hpp"
#include <limits>
#include <Eigen/Core>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#pragma GCC diagnostic ignored "-Wdeprecated-copy" // hurts, but is necessary
#include <pybind11/eigen.h>
#pragma GCC diagnostic pop

namespace py = pybind11;


namespace _pyves
{
    struct Particle
    {
        CARTESIAN position;
        CARTESIAN orientation;
        REAL sigma = std::numeric_limits<REAL>::signaling_NaN();
        REAL epsilon = std::numeric_limits<REAL>::signaling_NaN();
        REAL kappa = std::numeric_limits<REAL>::signaling_NaN();
        REAL gamma = std::numeric_limits<REAL>::signaling_NaN();
        std::string name = "UNDEFINED";
        
        Particle() = default;
        Particle(const Particle&) = default;
        // Particle(Particle&&) = default;
        Particle(CARTESIAN_CREF _pos, CARTESIAN_CREF _o) : position(_pos) , orientation(_o.normalized()) {}

        inline const auto getx() const { return position(0); }
        inline const auto gety() const { return position(1); }
        inline const auto getz() const { return position(2); }

        inline void setx(CARTESIAN::Scalar s) { position(0) = s; }
        inline void sety(CARTESIAN::Scalar s) { position(1) = s; }
        inline void setz(CARTESIAN::Scalar s) { position(2) = s; }

        inline const auto getux() const { return orientation(0); }
        inline const auto getuy() const { return orientation(1); }
        inline const auto getuz() const { return orientation(2); }

        inline auto getOrientation() const { return CARTESIAN_CREF(orientation); }
        inline void setOrientation(CARTESIAN_CREF v) { orientation = v; orientation.normalize(); }

        inline bool operator==(const Particle& other) const { return std::addressof(*this) == std::addressof(other); };
        inline bool operator!=(const Particle& other) const { return std::addressof(*this) != std::addressof(other); };
        
        inline std::string repr() const
        {
            return 
                std::string("<Particle ") + name + " (REAL=" + type_name<REAL>() + 
                ") at (" + 
                std::to_string(getx()) + "|" + std::to_string(gety()) + "|" + std::to_string(getz()) + 
                ") in (" + 
                std::to_string(getux()) + "|" + std::to_string(getuy()) + "|" + std::to_string(getuz()) + 
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
            .def_property("orientation", &Particle::getOrientation, &Particle::setOrientation)
            .def_readwrite("sigma", &Particle::sigma)
            .def_readwrite("epsilon", &Particle::epsilon)
            .def_readwrite("kappa", &Particle::kappa)
            .def_readwrite("gamma", &Particle::gamma)
            .def_readwrite("name", &Particle::name)
            .def_property("x", &Particle::getx, &Particle::setx)//, py::return_value_policy::reference_internal)
            .def_property("y", &Particle::gety, &Particle::sety)//, py::return_value_policy::reference_internal)
            .def_property("z", &Particle::getz, &Particle::setz)//, py::return_value_policy::reference_internal)
            .def_property_readonly("ux", &Particle::getux)//, &Particle::setux)//, py::return_value_policy::reference_internal)
            .def_property_readonly("uy", &Particle::getuy)//, &Particle::setuy)//, py::return_value_policy::reference_internal)
            .def_property_readonly("uz", &Particle::getuz)//, &Particle::setuz)//, py::return_value_policy::reference_internal)


            .def("__repr__", &Particle::repr);
            ;
    }
}