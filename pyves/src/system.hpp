#pragma once
#include "box.hpp"
// #include "parameters.hpp"
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
PYBIND11_MAKE_OPAQUE(std::vector<_pyves::Cell>)

namespace _pyves
{
    typedef std::vector<_pyves::Particle> ParticleContainer;
    typedef std::vector<_pyves::Cell> CellContainer;
    
    struct System
    {
        Box<PBC::ON> box;
        ParticleContainer particles;
        CellContainer cells;
    
        std::size_t cores = make_nan<std::size_t>();
        std::size_t threads = make_nan<std::size_t>();

        REAL temperature = make_nan<REAL>();
        REAL translation_min = make_nan<REAL>();
        REAL translation_max = make_nan<REAL>();
        REAL rotation_min = make_nan<REAL>();
        REAL rotation_max = make_nan<REAL>();
        std::size_t time_max = make_nan<std::size_t>();

        bool particleIsFree(const Particle&, REAL cutoff);
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



	    py::bind_vector<CellContainer>(m, "CellContainer" )
            .def(py::init<>())
            // .def("clear", &ParticleContainer::clear)
            // .def("pop_back", &ParticleContainer::pop_back)
            .def("__len__", [](const CellContainer& v) { return v.size(); })
            .def("__iter__", [](CellContainer& v) 
            {
                return py::make_iterator(std::begin(v), std::end(v));
            }, py::keep_alive<0, 1>())
            .def("__repr__", [](const CellContainer& v) {
                return "CellContainer\n[\n"
                    + std::accumulate(std::begin(v), std::end(v), std::string(""), [](std::string s, const Cell& p) 
                    { 
                        return s+"\t"+p.repr()+",\n";
                    })
                    + "]";
            })
            ;



        py::class_<System>(m, "System", py::dynamic_attr())
            .def(py::init<>())
            .def_readwrite("cores", &System::cores, py::return_value_policy::reference_internal)
            .def_readwrite("threads", &System::threads, py::return_value_policy::reference_internal)
            .def_readwrite("temperature", &System::temperature, py::return_value_policy::reference_internal)
            .def_readwrite("translation_min", &System::translation_min, py::return_value_policy::reference_internal)
            .def_readwrite("translation_max", &System::translation_max, py::return_value_policy::reference_internal)
            .def_readwrite("rotation_min", &System::rotation_min, py::return_value_policy::reference_internal)
            .def_readwrite("rotation_max", &System::rotation_max, py::return_value_policy::reference_internal)
            .def_readwrite("box", &System::box, py::return_value_policy::reference_internal)
            .def_readwrite("particles", &System::particles, py::return_value_policy::reference_internal)
            .def_readwrite("cells", &System::cells, py::return_value_policy::reference_internal)
            .def("particleIsFree", &System::particleIsFree)
            .def("assertIntegrity", &System::assertIntegrity)
            ;
    }


    
    bool System::particleIsFree(const Particle& subject, REAL cutoff)
    {
        const REAL cutoff_squared = cutoff*cutoff;
        return std::all_of(std::begin(particles), std::end(particles), [&](const auto& p) 
        { 
            return subject == p ? true : box.squaredDistanceParticle(subject, p) > cutoff_squared; 
        });
    }



    bool System::assertIntegrity() const
    {
        return all(
            std::all_of(std::begin(particles), std::end(particles), [](const auto& p) { return p.assertIntegrity(); }),
            std::all_of(std::begin(cells), std::end(cells), [](const auto& c) { return c.assertIntegrity(); }),
            std::isfinite(cores),
            std::isfinite(threads),
            std::isfinite(temperature),
            std::isfinite(translation_min),
            std::isfinite(translation_max),
            std::isfinite(rotation_min),
            std::isfinite(rotation_max),
            std::isfinite(time_max)
        );
    }
}