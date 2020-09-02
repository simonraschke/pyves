#pragma once

#include "definitions.hpp"
#include "particle.hpp"
#include "random.hpp"
#include <memory>
#include <Eigen/Geometry>

#include <pybind11/pybind11.h>
namespace py = pybind11;

#pragma GCC diagnostic ignored "-Wdeprecated-declarations" // hurts, but is necessary
#include <pybind11/eigen.h>
#pragma GCC diagnostic pop



namespace _pyves 
{
    enum class PBC : bool 
    { 
        ON=true, 
        OFF=false}
    ;

    template<PBC P> class Box;

    struct Particle;
}



// a simulation box implementation
// represents a box from 0 to x,x,z
// is ParameterDependentComponent
//
// may eiter be PBC::ON to calculate with PBC boundary conditions
// or may be PBC::OFF to caclulate without PBC boundary conditions
namespace _pyves
{
    template<PBC P>
    class Box
    {
    public:
        Box() = default;
        Box(REAL, REAL, REAL);
        Box(CARTESIAN_CREF);

        // set x by mutableAccess to Parameter base class
        void setLengthX(REAL);

        // set y by mutableAccess to Parameter base class
        void setLengthY(REAL);

        // set z by mutableAccess to Parameter base class
        void setLengthZ(REAL);

        REAL getLengthX() const;
        REAL getLengthY() const;
        REAL getLengthZ() const;
        REAL getVolume() const;

        CARTESIAN getCenter() const;

        // calcaulates the distance vector of two Particles
        // depending on PBC ON or OFF
        // called from anywhere else
        CARTESIAN distanceVector(CARTESIAN_CREF, CARTESIAN_CREF) const;
        // implementation for Particle base class. calls CARTESIAN version
        CARTESIAN distanceVectorParticle(const Particle&, const Particle&) const;


        // distance
        // calls squaredDistance and calculates std::sqrt
        REAL distance(CARTESIAN_CREF, CARTESIAN_CREF) const;
        // implementation for Particle base class. calls CARTESIAN version
        REAL distanceParticle(const Particle&, const Particle&) const;


        // squared distance
        REAL squaredDistance(CARTESIAN_CREF, CARTESIAN_CREF) const;
        // implementation for Particle base class. calls CARTESIAN version
        REAL squaredDistanceParticle(const Particle&, const Particle&) const;


        // scales down any give coordinates into the simulation box
        CARTESIAN scaleToBox(CARTESIAN) const ;
        // implementation for Particle base class. calls CARTESIAN version
        CARTESIAN scaleToBoxParticle(const Particle&) const;


        // scales down any give coordinates into the simulation box
        // VMD simulation box is from -x/2 to x/2
        // scales accordingly
        CARTESIAN scaleToBoxForVMD(CARTESIAN) const ;
        // implementation for Particle base class. calls CARTESIAN version
        CARTESIAN scaleToBoxForVMD(const Particle&) const;

        // checks if bounding_box contains coordinates
        // if PBC::ON calls scaleToBox before
        bool contains(CARTESIAN_CREF) const;
        // implementation for Particle base class. calls CARTESIAN version
        bool containsParticle(const Particle&) const;


        // destroy if derived is destroyed
        virtual ~Box() = default;

        // check if all parameters are set to make bounding_box
        // necessary for contains(const Particle&)
        void check_for_aligned_box_setup();

        CARTESIAN randomPointInside() const;
    
    protected:
        REAL x {0};
        REAL y {0};
        REAL z {0};

    private:
        std::unique_ptr<Eigen::AlignedBox<REAL,3>> bounding_box {nullptr};
    };



    template<PBC P>
    Box<P>::Box(REAL _x, REAL _y, REAL _z)
        : x(_x)
        , y(_y)
        , z(_z)
    {
        bounding_box.reset(nullptr);
        CARTESIAN vec (x,y,z);
        bounding_box = std::make_unique<Eigen::AlignedBox<REAL,3>>( CARTESIAN::Zero(), vec );
    }



    template<PBC P>
    Box<P>::Box(CARTESIAN_CREF _c)
        : x(_c[0])
        , y(_c[1])
        , z(_c[2])
    {
        bounding_box.reset(nullptr);
        bounding_box = std::make_unique<Eigen::AlignedBox<REAL,3>>( CARTESIAN::Zero(), _c );
    }



    template<PBC P>
    void Box<P>::setLengthX(REAL l)
    {
        x = l;
        check_for_aligned_box_setup();
    }



    template<PBC P>
    void Box<P>::setLengthY(REAL l)
    {
        y = l;
        check_for_aligned_box_setup();
    }



    template<PBC P>
    void Box<P>::setLengthZ(REAL l)
    {
        z = l;
        check_for_aligned_box_setup();
    }



    template<PBC P>
    REAL Box<P>::getLengthX() const
    {
        return x;
    }



    template<PBC P>
    REAL Box<P>::getLengthY() const
    {
        return y;
    }



    template<PBC P>
    REAL Box<P>::getLengthZ() const
    {
        return z;
    }



    template<PBC P>
    REAL Box<P>::getVolume() const
    {
        return bounding_box->volume();
    }



    template<PBC P>
    CARTESIAN Box<P>::getCenter() const
    {
        return bounding_box->center();
    }



    template<PBC P>
    void Box<P>::check_for_aligned_box_setup()
    {
        bounding_box.reset(nullptr);
        CARTESIAN vec (x,y,z);
        bounding_box = std::make_unique<Eigen::AlignedBox<REAL,3>>( CARTESIAN::Zero(), vec );
    }



    template<>
    inline CARTESIAN Box<PBC::ON>::distanceVector(CARTESIAN_CREF c1, CARTESIAN_CREF c2) const
    {
        CARTESIAN distance_CARTESIAN = c2-c1;
        distance_CARTESIAN(0) = distance_CARTESIAN(0) - x * std::round(static_cast<REAL>(distance_CARTESIAN(0)/(x)));
        distance_CARTESIAN(1) = distance_CARTESIAN(1) - y * std::round(static_cast<REAL>(distance_CARTESIAN(1)/(y)));
        distance_CARTESIAN(2) = distance_CARTESIAN(2) - z * std::round(static_cast<REAL>(distance_CARTESIAN(2)/(z)));
        return distance_CARTESIAN;
    }



    template<>
    inline CARTESIAN Box<PBC::OFF>::distanceVector(CARTESIAN_CREF c1, CARTESIAN_CREF c2) const
    {
        return (c2-c1);
    }



    template<PBC P>
    inline CARTESIAN Box<P>::distanceVectorParticle(const Particle& p1, const Particle& p2) const
    {
        return distanceVector(p1.position,p2.position);
    }



    template<PBC P>
    inline REAL Box<P>::squaredDistance(CARTESIAN_CREF c1, CARTESIAN_CREF c2) const 
    {
        return distanceVector(c1,c2).squaredNorm();
    }



    template<PBC P>
    inline REAL Box<P>::squaredDistanceParticle(const Particle& p1, const Particle& p2) const 
    {
        return squaredDistance(p1.position,p2.position);
    }



    template<PBC P>
    inline REAL Box<P>::distance(CARTESIAN_CREF c1, CARTESIAN_CREF c2) const 
    {
        return std::sqrt(squaredDistance(c1,c2));
    }



    template<PBC P>
    inline REAL Box<P>::distanceParticle(const Particle& p1, const Particle& p2) const 
    {
        return distance(p1.position,p2.position);
    }



    template<PBC P>
    CARTESIAN Box<P>::scaleToBox(CARTESIAN c) const 
    {
        while( c(0) > x ) c(0) -= x;
        while( c(1) > y ) c(1) -= y;
        while( c(2) > z ) c(2) -= z;

        while( c(0) < 0.f ) c(0) += x;
        while( c(1) < 0.f ) c(1) += y;
        while( c(2) < 0.f ) c(2) += z;
        return c;
    }



    template<PBC P>
    CARTESIAN Box<P>::scaleToBoxParticle(const Particle& p) const 
    {
        return scaleToBox(p.position);
    }



    template<PBC P>
    CARTESIAN Box<P>::scaleToBoxForVMD(CARTESIAN c) const 
    {
        c(0) = c(0) - x * std::round(c(0)/(x));
        c(1) = c(1) - y * std::round(c(1)/(y));
        c(2) = c(2) - z * std::round(c(2)/(z));
        return c;
    }



    template<PBC P>
    CARTESIAN Box<P>::scaleToBoxForVMD(const Particle& p) const 
    {
        return scaleToBoxForVMD(p.position);
    }



    template<>
    inline bool Box<PBC::ON>::contains(CARTESIAN_CREF c) const 
    {
        assert(bounding_box);
        return bounding_box->contains(scaleToBox(c));
    }



    template<>
    inline bool Box<PBC::OFF>::contains(CARTESIAN_CREF c) const 
    {
        assert(bounding_box);
        return bounding_box->contains(c);
    }



    template<PBC P>
    inline bool Box<P>::containsParticle(const Particle& p) const 
    {
        return contains(p.position);
    }



    template<PBC P>
    inline CARTESIAN Box<P>::randomPointInside() const 
    {
        return CARTESIAN
        (
            _pyves::random<CARTESIAN::Scalar>()(0.f,getLengthX()),
            _pyves::random<CARTESIAN::Scalar>()(0.f,getLengthY()),
            _pyves::random<CARTESIAN::Scalar>()(0.f,getLengthZ())
        );
    }



    template<PBC P>
    void declare_box(py::module& m, std::string typestr)
    {
        // using rvp = py::return_value_policy;
        using Class = Box<P>;
        std::string pyclass_name = std::string("Box")+typestr;
        py::class_<Class>(m, pyclass_name.c_str())
            .def(py::init<>())
            .def(py::init<REAL,REAL,REAL>())
            .def(py::init<CARTESIAN_CREF>())
            .def_property("x", &Class::getLengthX, &Class::setLengthX)
            .def_property("y", &Class::getLengthY, &Class::setLengthY)
            .def_property("z", &Class::getLengthZ, &Class::setLengthZ)
            .def_property_readonly("center", &Class::getCenter)
            .def_property_readonly("volume", &Class::getVolume)
            .def("distanceVector", &Class::distanceVector)
            .def("distanceVectorParticle", &Class::distanceVectorParticle)
            .def("distance", &Class::distance)
            .def("distanceParticle", &Class::distanceParticle)
            .def("squaredDistance", &Class::squaredDistance)
            .def("squaredDistanceParticle", &Class::squaredDistanceParticle)
            .def("scaleToBox", &Class::scaleToBox)
            .def("scaleToBoxParticle", &Class::scaleToBoxParticle)
            .def("contains", &Class::contains)
            .def("containsParticle", &Class::containsParticle)
            .def("randomPointInside", &Class::randomPointInside)
            ;
    }

    
    inline void bind_box(py::module& m)
    {
        declare_box<PBC::ON>(m, "PBC");
        declare_box<PBC::OFF>(m, "NoPBC");
    }
} // namespace _pyves