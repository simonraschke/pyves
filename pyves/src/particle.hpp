#pragma once

#include "definitions.hpp"
#include "utility.hpp"
#include <limits>
#include <exception>
#include <cmath>
#include <numeric>
#include <Eigen/Core>
#if __cplusplus >= 202000L
    #include <format>
#endif

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#if defined(__clang__)
    #pragma clang diagnostic ignored "-Wdeprecated-copy"
    #include <pybind11/eigen.h>
    // #pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
    #pragma GCC diagnostic ignored "-Wdeprecated-copy"
    #include <pybind11/eigen.h>
    #pragma GCC diagnostic pop
#elif defined(_MSC_VER)
    #include <pybind11/eigen.h>
#endif

namespace py = pybind11;

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>
#include <memory>



namespace _pyves
{
    struct Particle;   
    typedef std::vector<std::reference_wrapper<Particle>> ParticleRefContainer;
}



PYBIND11_MAKE_OPAQUE(_pyves::ParticleRefContainer)



namespace _pyves
{
        
    struct Particle
    {
    // private:
        CARTESIAN position;
        CARTESIAN orientation;
        std::unique_ptr<CARTESIAN> initial_position = nullptr;
        std::unique_ptr<CARTESIAN> initial_orientation = nullptr;
        REAL position_bound_radius_squared = std::numeric_limits<REAL>::max();
        REAL orientation_bound_radiant = std::numeric_limits<REAL>::max();

    // public:
        REAL sigma = std::numeric_limits<REAL>::signaling_NaN();
        REAL epsilon = std::numeric_limits<REAL>::signaling_NaN();
        REAL kappa = std::numeric_limits<REAL>::signaling_NaN();
        REAL gamma = std::numeric_limits<REAL>::signaling_NaN();
        std::string name = "UNDEFINED";
        
        Particle() = default;
        Particle(const Particle&);
        Particle& operator=(const Particle&);
        Particle(CARTESIAN_CREF _pos, CARTESIAN_CREF _o);
        Particle(CARTESIAN_CREF pos, CARTESIAN_CREF orien, REAL sigma, REAL eps, REAL kappa, REAL gamma, std::string name);



        inline CARTESIAN::Scalar getx() const { return position(0); }
        inline CARTESIAN::Scalar gety() const { return position(1); }
        inline CARTESIAN::Scalar getz() const { return position(2); }
        
        inline auto getPosition() const { return CARTESIAN_CREF(position); }
        inline CARTESIAN getInitialPosition() const
        {
            return *initial_position;
        }

        inline bool trySetPosition(CARTESIAN_CREF v)
        {
            if(!initial_orientation)
            {
                initial_orientation = std::make_unique<CARTESIAN>(v.normalized());
            }
            if((v - (*initial_position)).squaredNorm() < position_bound_radius_squared)
            {
                position = v;
                return true;
            }
            else
            {
                return false;
            }
            
        }

        inline void setInitialPosition(CARTESIAN_CREF p)
        {
            initial_position = std::make_unique<CARTESIAN>(p);
        }

        inline CARTESIAN::Scalar getux() const { return orientation(0); }
        inline CARTESIAN::Scalar getuy() const { return orientation(1); }
        inline CARTESIAN::Scalar getuz() const { return orientation(2); }



        inline auto getOrientation() const { return CARTESIAN_CREF(orientation); }
        inline void setOrientation(CARTESIAN_CREF v) 
        {
            orientation = v; 
            orientation.normalize(); 
        }

        inline bool trySetOrientation(CARTESIAN_CREF v)
        {
            if(!initial_orientation)
            {
                initial_orientation = std::make_unique<CARTESIAN>(v.normalized());
            }
            if(std::acos(v.normalized().dot(*initial_orientation)) < orientation_bound_radiant)
            {
                setOrientation(v);
                return true;
            }
            else
            {
                return false;
            }
            
        }

        inline CARTESIAN getInitialOrientation() const
        {
            return *initial_orientation;
        }

        inline void setInitialOrientation(CARTESIAN_CREF o)
        {
            initial_orientation = std::make_unique<CARTESIAN>(o);
        }



        inline bool operator==(const Particle& other) const { return this == &other; };
        inline bool operator!=(const Particle& other) const { return !(this == &other); };

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
            std::stringstream ss;
            ss  << "<Particle " << name << " (REAL=" << type_name<REAL>() << ") at " << position.format(VECTORFORMAT)
                << " in " << orientation.format(VECTORFORMAT) << " direction>";
            return ss.str();
        }  
        
        inline std::string detailed_repr() const
        {
            return repr() + 
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
            .def("size", &ParticleRefContainer::size)
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
            // .def(py::init<Particle>())
            .def(py::init<CARTESIAN_CREF, CARTESIAN_CREF>())
            .def(py::init<CARTESIAN_CREF, CARTESIAN_CREF, REAL, REAL, REAL, REAL, std::string>(), 
                 py::arg("pos"), py::arg("orien"), py::arg("sigma"), py::arg("eps"), py::arg("kappa"), py::arg("gamma"), py::arg("name"))
            .def(py::self == py::self)
            .def(py::self != py::self)
            .def("assertIntegrity", &Particle::assertIntegrity)
            .def_property("position", &Particle::getPosition, &Particle::trySetPosition)
            .def_property("orientation", &Particle::getOrientation, &Particle::trySetOrientation)
            .def_property("initial_position", &Particle::getInitialPosition, &Particle::setInitialPosition)
            .def_property("initial_orientation", &Particle::getInitialOrientation, &Particle::setInitialOrientation)
            .def_readwrite("sigma", &Particle::sigma)
            .def_readwrite("epsilon", &Particle::epsilon)
            .def_readwrite("kappa", &Particle::kappa)
            .def_readwrite("gamma", &Particle::gamma)
            .def_readwrite("name", &Particle::name)
            .def_property_readonly("x", &Particle::getx)//, &Particle::setx)
            .def_property_readonly("y", &Particle::gety)//, &Particle::sety)
            .def_property_readonly("z", &Particle::getz)//, &Particle::setz)
            .def_property_readonly("ux", &Particle::getux)
            .def_property_readonly("uy", &Particle::getuy)
            .def_property_readonly("uz", &Particle::getuz)
            .def("__repr__", &Particle::repr);
            ;
    }



    inline Particle::Particle(const Particle& other)
        : position(other.position)
        , orientation(other.orientation)
        , position_bound_radius_squared(other.position_bound_radius_squared)
        , orientation_bound_radiant(other.orientation_bound_radiant)
        , sigma(other.sigma)
        , epsilon(other.epsilon)
        , kappa(other.kappa)
        , gamma(other.gamma)
        , name(other.name)
    {
        if (!other.initial_position)
        {
            CARTESIAN copy = CARTESIAN(*other.initial_position.get());
            initial_position = std::make_unique<CARTESIAN>(CARTESIAN(copy));
        }
        if (!other.initial_orientation)
        {
            CARTESIAN copy = CARTESIAN(*other.initial_orientation.get());
            initial_orientation = std::make_unique<CARTESIAN>(CARTESIAN(copy));
        }
    }



    inline Particle& Particle::operator=(const Particle& other)
    {
        position = other.position;
        orientation = other.orientation;
        position_bound_radius_squared = other.position_bound_radius_squared;
        orientation_bound_radiant = other.orientation_bound_radiant;
        sigma = other.sigma;
        epsilon = other.epsilon;
        kappa = other.kappa;
        gamma = other.gamma;
        name = other.name;
        
        if (!other.initial_position)
        {
            CARTESIAN copy = CARTESIAN(*other.initial_position.get());
            initial_position = std::make_unique<CARTESIAN>(CARTESIAN(copy));
        }
        if (!other.initial_orientation)
        {
            CARTESIAN copy = CARTESIAN(*other.initial_orientation.get());
            initial_orientation = std::make_unique<CARTESIAN>(CARTESIAN(copy));
        }
        return *this;
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