#pragma once

#include <random>
#include <type_traits>



namespace _pyves
{
    // The static random Engine to be allocated once
    // if the application is shared memory parallelized
    //      this will not yield the same random numbers on every pull
    // std::random_device is to generate a real random number for the seed
    static struct RandomEngineInit
    {
        // construct and set seed via std::random_device
        explicit RandomEngineInit();

        // construct and set seed manually
        RandomEngineInit(const int);

        // get the seed
        int getSeed() const;

        // the pseudo random engine to be used publicly
        std::mt19937_64 pseudo_engine {};

    private:
        std::random_device true_engine {};
        const int seed;

    // call via enhance::RandomEngine.pseudo_engine
    } thread_local RandomEngine;



    // partially specialized template function to get random numbers
    // template<typename T>
    // T random(const T, const T);

    // partially specialized template function to get random numbers
    template<typename T>
    struct random
    {

        template<typename U=T, typename std::enable_if<std::is_integral<U>::value, std::nullptr_t>::type = nullptr>
        inline const T operator()(const T a, const T b) const
        {
            std::uniform_int_distribution<T> dist(a,b);
            return dist(RandomEngine.pseudo_engine);
        }

        template<typename U=T, typename std::enable_if<std::is_floating_point<U>::value, std::nullptr_t>::type = nullptr>
        inline const T operator()(const T a, const T b) const
        {
            std::uniform_real_distribution<T> dist(a,b);
            return dist(RandomEngine.pseudo_engine);
        }
    };
}