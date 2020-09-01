#pragma once

#include "definitions.hpp"
#include "type_name.hpp"
#include <limits>
#include <Eigen/Core>
#include <filesystem>

#include <pybind11/pybind11.h>
// #include <pybind11/operators.h>

// #pragma GCC diagnostic ignored "-Wdeprecated-declarations" // hurts, but is necessary
// #include <pybind11/eigen.h>
// #pragma GCC diagnostic pop

namespace py = pybind11;


namespace _pyves
{
    namespace __detail
    {
        template<typename T>
        struct NaN_t
        {
            constexpr auto operator() () const -> T { return std::numeric_limits<T>::signaling_NaN(); }
        };



        //
        // HARDWARE PARAMETER SET
        //
        struct HardwareParameters
        {
            std::size_t cores = NaN_t<std::size_t>()();
            std::size_t threads = NaN_t<std::size_t>()();
        };



        inline void bind_hardware_parameters(py::module& m)
        {                
            using Class = HardwareParameters;
            py::class_<Class>(m, "HardwareParameters", py::dynamic_attr())
                .def(py::init<>())
                .def(py::init<Class>())
                .def_readwrite("cores", &Class::cores)
                .def_readwrite("threads", &Class::threads)
                ;
        }



        //
        // IO PARAMETER SET
        //
        struct IOParameters
        {
            std::filesystem::path trajectory_path = ".";
            std::filesystem::path forcefield_path = ".";
            REAL time_delta = NaN_t<REAL>()();
        };



        inline void bind_io_parameters(py::module& m)
        {                
            using Class = IOParameters;
            py::class_<Class>(m, "IOParameters", py::dynamic_attr())
                .def(py::init<>())
                .def(py::init<Class>())
                .def_readwrite("trajectory_path", &Class::trajectory_path)
                .def_readwrite("forcefield_path", &Class::forcefield_path)
                ;
        }



        //
        // SYSTEM PARAMETER SET
        //
        struct SystemParameters
        {
            REAL temperature = NaN_t<REAL>()();
            REAL boxx = NaN_t<REAL>()();
            REAL boxy = NaN_t<REAL>()();
            REAL boxz = NaN_t<REAL>()();
            std::size_t time_max = NaN_t<std::size_t>()();
            REAL translation_min = NaN_t<REAL>()();
            REAL translation_max = NaN_t<REAL>()();
            REAL rotation_min = NaN_t<REAL>()();
            REAL rotation_max = NaN_t<REAL>()();
        };



        inline void bind_system_parameters(py::module& m)
        {                
            using Class = SystemParameters;
            py::class_<Class>(m, "SystemParameters", py::dynamic_attr())
                .def(py::init<>())
                .def(py::init<Class>())
                .def_readwrite("temperature", &Class::temperature)
                ;
        }
    }



    struct Parameters
    {
        __detail::HardwareParameters hardware;
        __detail::IOParameters io;
        __detail::SystemParameters system;
    };



    inline void bind_parameters(py::module& m)
    {
        __detail::bind_hardware_parameters(m);
        __detail::bind_io_parameters(m);
        __detail::bind_system_parameters(m);

        using Class = Parameters;
        py::class_<Class>(m, "Parameters", py::dynamic_attr())
            .def(py::init<>())
            .def(py::init<Class>())
            .def_readwrite("hardware", &Class::hardware)
            .def_readwrite("system", &Class::system)
            ;
    }
}