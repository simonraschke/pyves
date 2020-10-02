#pragma once

#include "definitions.hpp"
#include "box.hpp"
#include "math_utility.hpp"
#include <cmath>
#include <Eigen/Core>

#include <pybind11/pybind11.h>
namespace py = pybind11;


namespace _pyves
{
    struct Particle;
}



namespace _pyves::interaction
{
    REAL surface_potential(const Particle& p, const Box<PBC::ON>& box, REAL surface_width, REAL cutoff);
    REAL external_potential(const Particle& p, const Box<PBC::ON>& box, REAL surface_width, REAL cutoff);

    void bind_external_potential(py::module& m);
}

#include "particle.hpp"