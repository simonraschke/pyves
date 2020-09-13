#pragma once

#include "definitions.hpp"
#include <cmath>


namespace _pyves
{
    inline bool acceptByMetropolis(const REAL energy_difference, const REAL temperature)
    {
    #ifndef NDEBUG
        const REAL exr = std::exp(-energy_difference/temperature);
        const REAL ran = random<REAL>()(0.f,1.f);
        const bool acc = exr > ran;
        // vesDEBUG("energy difference: "<< energy_difference << "  temp: " << temperature << "  exp: " << exr <<"  random: " << ran << "  accepted: " << std::boolalpha << acc )
        assert( energy_difference < 0.f ? acc : true);
        return acc;
    #else
        return energy_difference < 0.f ? true : std::exp(-energy_difference/temperature) > random<REAL>(0.0,1.0);
    #endif
    }
} // namespace _pyves