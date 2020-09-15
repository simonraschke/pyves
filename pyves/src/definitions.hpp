#pragma once

#include <exception>

#if defined(__clang__)
    #pragma clang diagnostic ignored "-Wdeprecated-copy"
    #include <Eigen/Core>
    // #pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
    #pragma GCC diagnostic ignored "-Wdeprecated-copy"
    #include <Eigen/Core>
    #pragma GCC diagnostic pop
#elif defined(_MSC_VER)
    #include <Eigen/Core>
#endif

namespace _pyves
{
    typedef float REAL;
    typedef Eigen::Matrix<REAL, 3, 1> CARTESIAN;
    typedef Eigen::Ref<CARTESIAN> CARTESIAN_REF;
    typedef const Eigen::Ref<const CARTESIAN>& CARTESIAN_CREF;

    constexpr const double PI_4  = 0.785398163397448309616;
    constexpr const double PI_2  = 1.57079632679489661923;
    constexpr const double PI    = 3.14159265358979323846;
    constexpr const double TWOPI = 3.14159265358979323846*2;
}

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)


// #include <csignal>
// #include <iostream>
// #ifndef NDEBUG
//     #define vesDEBUG(x) {std::cerr << "[DEBUG] "; do { std::cerr << x; } while (0); std::cerr << '\n';}
// #else
//     #define vesDEBUG(x)
// #endif
// #define vesLOG(x) {std::clog << "[LOG] "; do { std::clog << x; } while (0); std::clog << '\n';}
// #define vesWARNING(x) {std::clog << "[WARNING] "; do { std::clog << x; } while (0); std::clog << '\n';}
// #define vesCRITICAL(x) {std::cerr << "[ERROR] "<< __FILE__ <<":" << __LINE__ << "  "; do { std::cerr << x; } while (0); std::cerr <<" raising SIGABRT\n"; std::exit(SIGABRT);}


//-------------------Eigen::IOFormat( prec, flag,                 coeffSep, rowSep, rowPre, rowSuf, matPre, matSuf )
#define ROWFORMAT    Eigen::IOFormat( 3,    Eigen::DontAlignCols, ", ",     " ",    " ",    "",     " ",    " " )
#define PYTHONFORMAT Eigen::IOFormat( 4,    0,                    ", ",     "\n",   "[",    "]",    "[",    "]" )
#define VECTORFORMAT Eigen::IOFormat( 3,    0,                    "|",      "|",    "",     "",     "(",    ")" )



#include <iostream>
template<typename T>
auto operator<<(std::ostream& os, const T& t) -> decltype(t.print(os), os) 
{ 
    t.print(os);
    return os;
}