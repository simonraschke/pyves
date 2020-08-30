#include "particle.hpp"
#include "system.hpp"

PYBIND11_MODULE(_pyves, m)
{
    _pyves::bind_particle(m);
    _pyves::bind_box(m);
}
