#pragma once

#include <random>
#include <type_traits>
#include <vector>
#include <deque>
#include <set>
#include <list>



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
    template<typename T, std::enable_if_t<std::is_integral<T>::value, bool> = true>
    inline T random(const T a, const T b)
    {
        std::uniform_int_distribution<T> dist(a,b);
        return dist(RandomEngine.pseudo_engine);
    }



    template<typename T, std::enable_if_t<std::is_floating_point<T>::value, bool> = true>
    inline T random(const T a, const T b)
    {
        std::uniform_real_distribution<T> dist(a,b);
        return dist(RandomEngine.pseudo_engine);
    }



    template<typename T, typename = void>
    struct is_iterator
    {
        static constexpr bool value = false;
    };

    template<typename T>
    struct is_iterator<T, std::enable_if_t<!std::is_same<typename std::iterator_traits<T>::value_type, void>::value>>
    {
        static constexpr bool value = true;
    };

    template<typename Iter, std::enable_if_t<is_iterator<Iter>::value, bool> = true>
    inline Iter random(Iter start, Iter end)
    {
        std::uniform_int_distribution<std::size_t> dist(0, std::distance(start, end) - 1);
        std::advance(start, dist(RandomEngine.pseudo_engine));
        return start;
    }



    template <typename T>
    struct is_container 
    {
        static constexpr bool value = false;
    };

    template <typename T,typename Alloc>
    struct is_container<std::vector<T,Alloc> > 
    {
        static const bool value = true;
    };

    template <typename T,typename Alloc>
    struct is_container<std::deque<T,Alloc> > 
    {
        static const bool value = true;
    };

    template <typename T,typename Alloc>
    struct is_container<std::set<T,Alloc> > 
    {
        static const bool value = true;
    };

    template <typename T,typename Alloc>
    struct is_container<std::list<T,Alloc> > 
    {
        static const bool value = true;
    };

    template <typename Container, std::enable_if_t<is_container<Container>::value, bool> = true>
    inline auto random(Container& c) -> decltype(*begin(c))& 
    {
        return *random(std::begin(c), std::end(c));
    }
}