#pragma once

#include "type_name.hpp"
#include "point.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;



template<typename T>
struct Position : public Point<T>
{
    Position(): Point<T>()
    {
        ;
    }
    
    Position(Point<T> _p): Point<T>(_p)
    {
        ;
    }

    Position(T _x, T _y, T _z): Point<T>(_x,_y,_z)
    {
        ;
    }

    // using Point<T>::Point;

    inline void translation(const Point<T>& other)
    {
        this->operator+=(other);
    }

    inline constexpr std::string repr() const
    {
        return std::string("<Position<") + type_name<T>() + "> at " + std::to_string(this->x) + " " + std::to_string(this->y) + " " + std::to_string(this->z) + ">" ;
    }
};



template<typename T>
void declare_position(py::module &m, std::string typestr)
{
    using Class = Position<T>;
    std::string pyclass_name = std::string("Position") + typestr;
    py::class_<Class, Point<T>>(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
        .def(py::init<>())
        .def(py::init<Point<T>>())
        .def(py::init<T,T,T>())
        .def("translation", &Class::translation)
        .def("__repr__", &Class::repr);
}



inline void bind_position(py::module &m) 
{
    declare_position<float>(m, "f");
    declare_position<double>(m, "d");
    declare_position<int>(m, "i");
    declare_position<long>(m, "l");
}

// /// Arbitrary axis angle rotation
// Matrix Matrix::Rotation( const Vector & u, real radians )
// {
//   real c = cos( radians ) ;
//   real l_c = 1 - c ;

//   real s = sin( radians ) ;
  
//   //ROW MAJOR
//   return Matrix(
//     u.x*u.x + (1 - u.x*u.x)*c,      u.x*u.y*l_c + u.z*s,        u.x*u.z*l_c - u.y*s,  0,
//           u.x*u.y*l_c - u.z*s,  u.y*u.y+(1 - u.y*u.y)*c,        u.y*u.z*l_c + u.x*s,  0,
//           u.x*u.z*l_c + u.y*s,      u.y*u.z*l_c - u.x*s,  u.z*u.z + (1 - u.z*u.z)*c,  0,
//                             0,                        0,                          0,  1
//   ) ;
// }