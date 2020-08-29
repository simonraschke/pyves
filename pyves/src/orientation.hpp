#pragma once

#include "type_name.hpp"
#include "point.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;



template<typename T>
struct Orientation : public Point<T>
{
    Orientation(): Point<T>()
    {
        ;
    }

    Orientation(T _x, T _y, T _z): Point<T>(_x,_y,_z)
    {
        ;
    }

    inline Orientation<T> operator+(const Point<T>& other) const
    {
        return Orientation<T>(this->x+other.x, this->y+other.y, this->z+other.z);
    }
    
    inline void translation(const Point<T>& other)
    {
        this->operator+=(other);
    }

    inline constexpr std::string repr() const
    {
        return std::string("<Orientation<") + type_name<T>() + "> at " + std::to_string(this->x) + " " + std::to_string(this->y) + " " + std::to_string(this->z) + ">" ;
    }
};



template<typename T>
void declare_orientation(py::module &m, std::string typestr)
{
    using Class = Orientation<T>;
    std::string pyclass_name = std::string("Orientation") + typestr;
    py::class_<Class, Point<T>>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
        .def(py::init<>())
        .def(py::init<T,T,T>())
        .def("translation", &Class::translation)
        .def("__repr__", &Class::repr);
}



inline void bind_orientation(py::module &m) 
{
    declare_orientation<float>(m, "f");
    declare_orientation<double>(m, "d");
    declare_orientation<int>(m, "i");
    declare_orientation<long>(m, "l");
}
