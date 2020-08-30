#include "catch.hpp"
#include "particle.hpp"



using namespace _pyves;



TEST_CASE("Particle Constructor Test")
{
    Particle p1;
    CHECK_NOTHROW(Particle());
    CHECK_NOTHROW(Particle(p1));
    CHECK_NOTHROW(Particle(CARTESIAN(1,1,1), CARTESIAN(1,0,0)));
}



TEST_CASE("Particle Member Test")
{
    SECTION("position and orientation access")
    {
        Particle p1 = Particle(CARTESIAN(4, 3.2, -4.5), CARTESIAN(1,0,0));
        CHECK(p1.position(0) == Approx(static_cast<REAL>(   4)));
        CHECK(p1.position(1) == Approx(static_cast<REAL>( 3.2)));
        CHECK(p1.position(2) == Approx(static_cast<REAL>(-4.5)));
        CHECK(p1.orientation(0) == Approx(static_cast<REAL>(1)));
        CHECK(p1.orientation(1) == Approx(static_cast<REAL>(0)));
        CHECK(p1.orientation(2) == Approx(static_cast<REAL>(0)));

        p1.position = {1.2, -4, 3.1};
        CHECK(p1.position(0) == Approx(static_cast<REAL>(1.2)));
        CHECK(p1.position(1) == Approx(static_cast<REAL>( -4)));
        CHECK(p1.position(2) == Approx(static_cast<REAL>(3.1)));

        p1.orientation = {-1,0,1};
        CHECK(p1.orientation(0) == Approx(static_cast<REAL>(-1)));
        CHECK(p1.orientation(1) == Approx(static_cast<REAL>( 0)));
        CHECK(p1.orientation(2) == Approx(static_cast<REAL>( 1)));
    }


    SECTION("operators")
    {
        Particle p1(CARTESIAN(1,1,1), CARTESIAN(1,0,0));
        Particle p2(CARTESIAN(2,-2,0.1), CARTESIAN(1,0,0));
        CHECK(p1 == p1);
        CHECK(p1 != p2);
        CHECK_FALSE(p1 == p2);
        CHECK_FALSE(p1 != p1);
    }


    SECTION("simple additions")
    {
        Particle p1(CARTESIAN(1,1,1), CARTESIAN(1,0,0));
        Particle p2(CARTESIAN(2,-2,0.1), CARTESIAN(1,0,0));
        Particle p3;
        p3.position = p1.position + p2.position;
        CHECK(p3.position(0) == Approx(static_cast<REAL>(3)));
        CHECK(p3.position(1) == Approx(static_cast<REAL>(-1)));
        CHECK(p3.position(2) == Approx(static_cast<REAL>(1.1)));

        p3.position += p1.position;
        CHECK(p3.position(0) == Approx(static_cast<REAL>(4)));
        CHECK(p3.position(1) == Approx(static_cast<REAL>(0)));
        CHECK(p3.position(2) == Approx(static_cast<REAL>(2.1)));

        p3.position += CARTESIAN(1,1,1);
        CHECK(p3.position(0) == Approx(static_cast<REAL>(5)));
        CHECK(p3.position(1) == Approx(static_cast<REAL>(1)));
        CHECK(p3.position(2) == Approx(static_cast<REAL>(3.1)));
    }


    SECTION("direct x,y,z, ux,uy,uz access")
    {
        Particle p;
        CHECK(p.x() == Approx(static_cast<REAL>(0)).margin(1e-7));
        CHECK(p.y() == Approx(static_cast<REAL>(0)).margin(1e-7));
        CHECK(p.z() == Approx(static_cast<REAL>(0)).margin(1e-7));
        p.x() = 1;
        p.y() = 2.2;
        p.z() = -1.3;
        CHECK(p.x() == Approx(static_cast<REAL>(1)));
        CHECK(p.y() == Approx(static_cast<REAL>(2.2)));
        CHECK(p.z() == Approx(static_cast<REAL>(-1.3)));
        CHECK(p.position(0) == Approx(static_cast<REAL>(1)));
        CHECK(p.position(1) == Approx(static_cast<REAL>(2.2)));
        CHECK(p.position(2) == Approx(static_cast<REAL>(-1.3)));
    }
}
