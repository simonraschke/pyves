#pragma once

#include "definitions.hpp"
#include "utility.hpp"
#include "box.hpp"
#include "interaction.hpp"
#include "external_potential.hpp"
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
        CARTESIAN position;
        CARTESIAN orientation;
        std::unique_ptr<CARTESIAN> initial_position = nullptr;
        std::unique_ptr<CARTESIAN> initial_orientation = nullptr;

        REAL sigma = std::numeric_limits<REAL>::signaling_NaN();
        REAL epsilon = std::numeric_limits<REAL>::signaling_NaN();
        REAL kappa = std::numeric_limits<REAL>::signaling_NaN();
        REAL gamma = std::numeric_limits<REAL>::signaling_NaN();
        REAL position_bound_radius_squared = std::numeric_limits<REAL>::max();
        REAL orientation_bound_radiant = std::numeric_limits<REAL>::max();
        
        REAL surface_affinity_translation = 0;
        REAL surface_affinity_rotation = 0;

        REAL self_affinity = 1.0;
        REAL other_affinity = 1.0;

        std::string name = "UNDEFINED";

        ParticleRefContainer neighbors;
        
        Particle() = default;
        Particle(const Particle&);
        Particle(CARTESIAN_CREF pos, CARTESIAN_CREF orien);
        Particle(CARTESIAN_CREF pos, CARTESIAN_CREF orien, REAL sigma, REAL eps, REAL kappa, REAL gamma, std::string name);
        Particle& operator=(const Particle&);

        CARTESIAN::Scalar getx() const;
        CARTESIAN::Scalar gety() const;
        CARTESIAN::Scalar getz() const;
        CARTESIAN::Scalar getux() const;
        CARTESIAN::Scalar getuy() const;
        CARTESIAN::Scalar getuz() const;

        CARTESIAN getPosition() const;
        CARTESIAN getInitialPosition() const;

        bool trySetPosition(CARTESIAN_CREF);
        void setInitialPosition(CARTESIAN);
        void setPosition(CARTESIAN_CREF);

        CARTESIAN getOrientation() const;
        CARTESIAN getInitialOrientation() const;

        bool trySetOrientation(CARTESIAN_CREF);
        void setInitialOrientation(CARTESIAN);
        void setOrientation(CARTESIAN_CREF);

        void updateNeighborList(const ParticleRefContainer&, const Box<PBC::ON>&, REAL);

        REAL potentialEnergy(const Box<PBC::ON>&, REAL) const;
        REAL chi(const Box<PBC::ON>&, REAL) const;
        REAL externalPotential(const Box<PBC::ON>&, REAL cutoff, REAL surface_width) const;

        bool operator==(const Particle& other) const;
        bool operator!=(const Particle& other) const;

        bool assertIntegrity() const;
        
        std::string repr() const;
        std::string detailed_repr() const;
    };



    void bind_particle(py::module&);
}