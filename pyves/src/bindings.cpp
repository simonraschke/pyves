#include "point.hpp"
// #include "position.hpp"
#include "particle.hpp"
// #include "orientation.hpp"

PYBIND11_MODULE(_pyves, m)
{
    // m.doc();
    bind_point(m);
    // bind_position(m);
    // bind_orientation(m);
    _pyves::bind_particle(m);
}
