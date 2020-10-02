#include "particle.hpp"
#include "cell.hpp"
#include "stepwidth_alignment_unit.hpp"
#include "system.hpp"
#include "interaction.hpp"
#include "external_potential.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;

PYBIND11_MODULE(_pyves, m)
{
    _pyves::bind_particle(m);
    _pyves::bind_box(m);
    _pyves::bind_cell(m);
    _pyves::bind_sw_alignment_unit(m);
    _pyves::bind_system(m);
    _pyves::interaction::bind_interaction(m);
    _pyves::interaction::bind_external_potential(m);

    m.def("concurrency_model", []() -> py::str {
        #ifdef PYVES_USE_TBB
        return py::str("TBB");
        #else
        return py::str("Taskflow");
        #endif
    });
}
