#pragma once
#include "definitions.hpp"
#include "box.hpp"
#include "parameters.hpp"
#include "particle.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <pybind11/stl.h>
#include <pybind11/buffer_info.h>
#include <vector>

PYBIND11_MAKE_OPAQUE(std::vector<_pyves::Particle, std::allocator<_pyves::Particle>>);

using ParticleContainer = std::vector<_pyves::Particle, std::allocator<_pyves::Particle>>;

namespace _pyves
{
    struct System
    {
        Box<PBC::ON> box;
        Parameters prms;
        ParticleContainer particles;
        // cells
        // 
        

        bool assertIntegrity() const;
    };



    inline void bind_system(py::module& m)
    {
        py::class_<System>(m, "System", py::dynamic_attr())
            .def(py::init<>())
            .def_readwrite("prms", &System::prms)
            .def_readwrite("particles", &System::particles)
            ;
        // py::bind_vector<std::vector<Particle>>(m, "System.particles");
    }
}