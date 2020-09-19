#pragma once

#include "definitions.hpp"
#include "box.hpp"
#include <cmath>
#include <unordered_map>
#include <Eigen/Core>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
namespace py = pybind11;


namespace _pyves
{
    struct Particle;

    namespace __detail
    {
        template<class A, class B>
        using map_t = std::unordered_map<A,B>;
        
    }
    typedef __detail::map_t<std::string, __detail::map_t<std::string, __detail::map_t<std::string, REAL>>> LookupTable_t;
}



PYBIND11_MAKE_OPAQUE(_pyves::LookupTable_t)



namespace _pyves
{
    REAL interaction(const Particle& p1, const Particle& p2, const Box<PBC::ON>& box, REAL cutoff);
    REAL interactionWithLookup(const Particle& p1, const Particle& p2, const Box<PBC::ON>& box, REAL cutoff, const LookupTable_t& table);

    REAL calculateInteractionOptimumA(const Particle& p1, const Particle& p2);
    REAL calculateInteractionOptimumB(const Particle& p1, const Particle& p2);
    REAL calculateInteractionOptimumC1(const Particle& p1, const Particle& p2);
    REAL calculateInteractionOptimumC2(const Particle& p1, const Particle& p2);

    LookupTable_t calculateInteractionLookupTable(std::vector<Particle> particles);



    void bind_interaction(py::module& m);
}

#include "particle.hpp"