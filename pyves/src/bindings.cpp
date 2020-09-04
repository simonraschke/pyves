// #include "parameters.hpp"
#include "particle.hpp"
#include "cell.hpp"
#include "system.hpp"
#include "interaction.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;

PYBIND11_MODULE(_pyves, m)
{
    _pyves::bind_particle(m);
    _pyves::bind_box(m);
    _pyves::bind_cell(m);
    _pyves::bind_system(m);
    _pyves::bind_interaction(m);
}
