#include "random.hpp"


namespace _pyves
{
    RandomEngineInit::RandomEngineInit()
        : seed(true_engine())
    {
        pseudo_engine.seed(seed);
    }



    RandomEngineInit::RandomEngineInit(const int __seed)
        : seed(__seed)
    {
        pseudo_engine.seed(seed);
    }



    int RandomEngineInit::getSeed() const
    {
        return seed;
    }
}
