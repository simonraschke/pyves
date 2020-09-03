#pragma once
#include "box.hpp"
#include "parameters.hpp"
#include "particle.hpp"
#include "cell.hpp"
#include "utility.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>
#include <numeric>

PYBIND11_MAKE_OPAQUE(std::vector<_pyves::Particle>)

namespace _pyves
{
    typedef std::vector<_pyves::Particle> ParticleContainer;
    
    struct System
    {
        Box<PBC::ON> box;
        Parameters prms;
        ParticleContainer particles = {};
        // cells
        // 
        
        void prepare_simulation();
        bool assertIntegrity() const;
    };



    inline void bind_system(py::module& m)
    {
	    py::bind_vector<ParticleContainer>( m, "ParticleContainer" )
            .def(py::init<>())
            // .def("clear", &ParticleContainer::clear)
            // .def("pop_back", &ParticleContainer::pop_back)
            .def("__len__", [](const ParticleContainer& v) { return v.size(); })
            .def("__iter__", [](ParticleContainer& v) 
            {
                return py::make_iterator(std::begin(v), std::end(v));
            }, py::keep_alive<0, 1>())
            .def("__repr__", [](const ParticleContainer& v) {
                return "ParticleContainer\n[\n"
                    + std::accumulate(std::begin(v), std::end(v), std::string(""), [](std::string s, const Particle& p) 
                    { 
                        return s+"\t"+p.repr()+",\n";
                    })
                    + "]";
            })
            ;

        py::class_<System>(m, "System", py::dynamic_attr())
            .def(py::init<>())
            .def_readwrite("prms", &System::prms)
            .def_readwrite("particles", &System::particles)
            .def("assertIntegrity", &System::assertIntegrity)
            ;
    }



    bool System::assertIntegrity() const
    {
        return all(
            std::all_of(std::begin(particles), std::end(particles), [](const auto& p) { return p.assertIntegrity(); })
        );
    }
}