#pragma once

#include "type_name.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

namespace py = pybind11;



template<typename T>
struct Point{
    T x = {0};
    T y = {0};
    T z = {0};

    Point() {}
    Point(T _x, T _y, T _z) : x(_x), y(_y), z(_z)  {}
    virtual ~Point() = default;

    inline constexpr std::string repr() const
    {
        return std::string("<Point<") + type_name<T>() + "> at " + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + ">" ;
    }

    inline Point<T> operator+(const Point<T>& other) const
    {
        return Point<T>(x+other.x, y+other.y, z+other.z);
    }

    inline Point<T> operator-(const Point<T>& other) const
    {
        return Point<T>(x-other.x, y-other.y, z-other.z);
    }
    
    inline Point<T>& operator+=(const Point<T>& other)
    {
        this->x += other.x;
        this->y += other.y;
        this->z += other.z;
        return *this;
    }
    
    inline T dot(const Point<T>& other) const
    {
        return x*other.x + y*other.y + z*other.z;
    }
};



template<typename T>
void declare_point(py::module &m, std::string typestr)
 {
    using Class = Point<T>;
    std::string pyclass_name = std::string("Point") + typestr;
    py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
        .def(py::init<>())
        .def(py::init<T,T,T>())
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self += py::self)
        .def("dot", &Class::dot)
        .def_readwrite("x", &Class::x)
        .def_readwrite("y", &Class::y)
        .def_readwrite("z", &Class::z)
        .def("__repr__", &Class::repr);
}



inline void bind_point(py::module &m) 
{
    declare_point<float>(m, "f");
    declare_point<double>(m, "d");
    declare_point<int>(m, "i");
    declare_point<long>(m, "l");
}
