#include "definitions.hpp"
#include "math_utility.hpp"
// #include "parallel.hpp"
#include <Eigen/Core>
#include <tbb/task_arena.h>

#include <pybind11/pybind11.h>
namespace py = pybind11;


namespace _pyves
{

    struct Particle
    {
        Eigen::Array<REAL,3,1> coordinates;
        Eigen::Array<REAL,3,1> orientation;

        Particle()
        {
            tbb::task_arena limited(3);
            
            limited.execute([]{std::cout << 1;});
            limited.execute([]{std::cout << 2;});
            limited.execute([]{std::cout << 3;});

        }
    };

}



inline void bind_particle(py::module& m) 
{
     py::class_<_pyves::Particle>(m, "Particle")
        .def(py::init<>())
        // .def(py::init<T,T,T>())
        // .def(py::init<Point<int>>())
        // .def(py::init<Point<long>>())
        // .def(py::init<Point<float>>())
        // .def(py::init<Point<double>>())
        // .def(py::self + py::self)
        // .def(py::self - py::self)
        // .def(py::self += py::self)
        // .def("dot", &Class::dot)
        // .def_readwrite("x", &Class::x)
        // .def_readwrite("y", &Class::y)
        // .def_readwrite("z", &Class::z)
        // .def("__repr__", &Class::repr);
        ;
}