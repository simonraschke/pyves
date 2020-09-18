#pragma once

#include "definitions.hpp"
#include "box.hpp"
#include <cmath>
#include <unordered_map>
#include <Eigen/Core>

#include <pybind11/pybind11.h>
namespace py = pybind11;


namespace _pyves
{
    struct Particle;

    namespace __detail
    {
        template<class A, class B>
        using map_t = std::map<A,B>;
        
    }
    typedef __detail::map_t<std::string, __detail::map_t<std::string, __detail::map_t<std::string, REAL>>> LookupTable_t;



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