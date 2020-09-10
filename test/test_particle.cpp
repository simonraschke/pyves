#include "catch.hpp"
#include "particle.hpp"



using namespace _pyves;



TEST_CASE("Particle Constructor Test")
{
    Particle p1;
    CHECK_NOTHROW(Particle());
    // CHECK_NOTHROW(Particle(p1));
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
        CHECK(p1.getOrientation()(0) == Approx(static_cast<REAL>(1)));
        CHECK(p1.getOrientation()(1) == Approx(static_cast<REAL>(0)));
        CHECK(p1.getOrientation()(2) == Approx(static_cast<REAL>(0)));

        p1.position = CARTESIAN(1.2, -4, 3.1);
        CHECK(p1.position(0) == Approx(static_cast<REAL>(1.2)));
        CHECK(p1.position(1) == Approx(static_cast<REAL>( -4)));
        CHECK(p1.position(2) == Approx(static_cast<REAL>(3.1)));

        p1.setOrientation(CARTESIAN(-1,0,1));
        // p1.orientation = {-1,0,1};
        CHECK(p1.getOrientation()(0) == Approx(static_cast<REAL>(-1/std::sqrt(2))));
        CHECK(p1.getOrientation()(1) == Approx(static_cast<REAL>( 0)));
        CHECK(p1.getOrientation()(2) == Approx(static_cast<REAL>( 1/std::sqrt(2))));
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
        Particle p1(CARTESIAN(1,1,1), CARTESIAN(8,0,0));
        Particle p2(CARTESIAN(2,-2,0.1), CARTESIAN(2,0,0));
        Particle p3;
        p3.position = p1.position + p2.position;
        // p3.setPosition(p1.position + p2.position);
        CHECK(p3.position(0) == Approx(static_cast<REAL>(3)));
        CHECK(p3.position(1) == Approx(static_cast<REAL>(-1)));
        CHECK(p3.position(2) == Approx(static_cast<REAL>(1.1)));

        p3.position += p1.position;
        // p3.setPosition(p3.position + p1.position);
        CHECK(p3.position(0) == Approx(static_cast<REAL>(4)));
        CHECK(p3.position(1) == Approx(static_cast<REAL>(0)));
        CHECK(p3.position(2) == Approx(static_cast<REAL>(2.1)));

        p3.position += CARTESIAN(1,1,1);
        // p3.setPosition(p3.position + CARTESIAN(1,1,1));
        CHECK(p3.position(0) == Approx(static_cast<REAL>(5)));
        CHECK(p3.position(1) == Approx(static_cast<REAL>(1)));
        CHECK(p3.position(2) == Approx(static_cast<REAL>(3.1)));
        
        CHECK(p1.getux() == Approx(1));
        CHECK(p1.getuy() == Approx(0));
        CHECK(p1.getuz() == Approx(0));
        CHECK(p2.getux() == Approx(1));
        CHECK(p2.getuy() == Approx(0));
        CHECK(p2.getuz() == Approx(0));
    }



    SECTION("direct x,y,z, ux,uy,uz access")
    {
        Particle p;
        CHECK(p.getx() == Approx(static_cast<REAL>(0)).margin(1e-7));
        CHECK(p.gety() == Approx(static_cast<REAL>(0)).margin(1e-7));
        CHECK(p.getz() == Approx(static_cast<REAL>(0)).margin(1e-7));
        p.position = CARTESIAN(1, 2.2, -1.3);
        // p.setx(1);
        // p.sety(2.2);
        // p.setz(-1.3);
        CHECK(p.getx() == Approx(static_cast<REAL>(1)));
        CHECK(p.gety() == Approx(static_cast<REAL>(2.2)));
        CHECK(p.getz() == Approx(static_cast<REAL>(-1.3)));
        CHECK(p.position(0) == Approx(static_cast<REAL>(1)));
        CHECK(p.position(1) == Approx(static_cast<REAL>(2.2)));
        CHECK(p.position(2) == Approx(static_cast<REAL>(-1.3)));
        p.setOrientation(CARTESIAN(2.2345234,0,0));
        CHECK(p.getux() == Approx(1));
        CHECK(p.getuy() == Approx(0));
        CHECK(p.getuz() == Approx(0));
    }
}