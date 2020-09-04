#pragma once

#include "definitions.hpp"
#include "utility.hpp"
#include <limits>
#include <exception>
#include <cmath>
#include <Eigen/Core>
#if __cplusplus >= 202000L
    #include <format>
#endif

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#pragma GCC diagnostic ignored "-Wdeprecated-copy" // hurts, but is necessary
#include <pybind11/eigen.h>
#pragma GCC diagnostic pop

namespace py = pybind11;

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>

namespace _pyves
{
    struct Particle;   
    typedef std::vector<std::reference_wrapper<Particle>> ParticleRefContainer;
}



PYBIND11_MAKE_OPAQUE(std::vector<std::reference_wrapper<_pyves::Particle>>)



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
        Particle(CARTESIAN_CREF _pos, CARTESIAN_CREF _o);
        Particle(CARTESIAN_CREF pos, CARTESIAN_CREF orien, REAL sigma, REAL eps, REAL kappa, REAL gamma, std::string name);

        inline CARTESIAN::Scalar getx() const { return position(0); }
        inline CARTESIAN::Scalar gety() const { return position(1); }
        inline CARTESIAN::Scalar getz() const { return position(2); }

        inline void setx(CARTESIAN::Scalar s) { position(0) = s; }
        inline void sety(CARTESIAN::Scalar s) { position(1) = s; }
        inline void setz(CARTESIAN::Scalar s) { position(2) = s; }

        inline CARTESIAN::Scalar getux() const { return orientation(0); }
        inline CARTESIAN::Scalar getuy() const { return orientation(1); }
        inline CARTESIAN::Scalar getuz() const { return orientation(2); }

        inline auto getOrientation() const { return CARTESIAN_CREF(orientation); }
        inline void setOrientation(CARTESIAN_CREF v) { orientation = v; orientation.normalize(); }

        inline bool operator==(const Particle& other) const { return std::addressof(*this) == std::addressof(other); };
        inline bool operator!=(const Particle& other) const { return std::addressof(*this) != std::addressof(other); };

        inline bool assertIntegrity() const
        {
            return all(
                all(std::isfinite(position[0]), std::isfinite(position[1]), std::isfinite(position[2])),
                all(std::isfinite(orientation[0]), std::isfinite(orientation[1]), std::isfinite(orientation[2])),
                std::isfinite(sigma),
                std::isfinite(epsilon),
                std::isfinite(kappa),
                std::isfinite(gamma),
                name != "UNDEFINED"
            );
        }
        
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
        
        inline std::string detailed_repr() const
        {
            return 
                std::string("<Particle ") + name + " (REAL=" + type_name<REAL>() + 
                ") at (" + 
                std::to_string(getx()) + "|" + std::to_string(gety()) + "|" + std::to_string(getz()) + 
                ") in (" + 
                std::to_string(getux()) + "|" + std::to_string(getuy()) + "|" + std::to_string(getuz()) + 
            #if __cplusplus <= 201703L
                string_format(") direction with sigma=%f eps=%f kappa=%f gamma=%f >", sigma, epsilon, kappa, gamma);
            #else
                std::format(") direction with sigma={} eps={} kappa={} gamma={} >", sigma, epsilon, kappa, gamma);
            #endif
        }  
    };



    inline void bind_particle(py::module& m) 
    {
        
        py::bind_vector<ParticleRefContainer>( m, "ParticleRefContainer" )
            .def(py::init<>())
            .def("__len__", [](const ParticleRefContainer& v) { return v.size(); })
            .def("__iter__", [](ParticleRefContainer& v) 
            {
                return py::make_iterator(std::begin(v), std::end(v));
            }, py::keep_alive<0, 1>())
            .def("__repr__", [](const ParticleRefContainer& v) {
                return "ParticleRefContainer\n[\n"
                    + std::accumulate(std::begin(v), std::end(v), std::string(""), [](std::string s, const Particle& p) 
                    { 
                        return s+"\t"+p.repr()+",\n";
                    })
                    + "]";
            })
            ;
        
        py::class_<Particle>(m, "Particle", py::dynamic_attr())
            .def(py::init<>())
            .def(py::init<Particle>())
            .def(py::init<CARTESIAN_CREF, CARTESIAN_CREF>())
            .def(py::init<CARTESIAN_CREF, CARTESIAN_CREF, REAL, REAL, REAL, REAL, std::string>(), 
                 py::arg("pos"), py::arg("orien"), py::arg("sigma"), py::arg("eps"), py::arg("kappa"), py::arg("gamma"), py::arg("name"))
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("assertIntegrity", &Particle::assertIntegrity)
            .def_readwrite("position", &Particle::position)
            .def_property("orientation", &Particle::getOrientation, &Particle::setOrientation)
            .def_readwrite("sigma", &Particle::sigma)
            .def_readwrite("epsilon", &Particle::epsilon)
            .def_readwrite("kappa", &Particle::kappa)
            .def_readwrite("gamma", &Particle::gamma)
            .def_readwrite("name", &Particle::name)
            .def_property("x", &Particle::getx, &Particle::setx)
            .def_property("y", &Particle::gety, &Particle::sety)
            .def_property("z", &Particle::getz, &Particle::setz)
            .def_property_readonly("ux", &Particle::getux)
            .def_property_readonly("uy", &Particle::getuy)
            .def_property_readonly("uz", &Particle::getuz)
            .def("__repr__", &Particle::repr);
            ;
    }



    inline Particle::Particle(CARTESIAN_CREF pos, CARTESIAN_CREF orien)
        : position(pos)
        , orientation(orien.normalized()) 
    {

    }

    

    inline Particle::Particle(CARTESIAN_CREF pos, CARTESIAN_CREF orien, REAL sigma, REAL eps, REAL kappa, REAL gamma, std::string name)
        : position(pos)
        , orientation(orien.normalized())
        , sigma(sigma)
        , epsilon(eps)
        , kappa(kappa)
        , gamma(gamma)
        , name(name)
    {
        // if(!assertIntegrity())
        // {
        //     throw std::logic_error(std::string("Particle should be completely initialized, but is not: ")+detailed_repr());
        // }
    }
}