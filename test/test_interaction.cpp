#include "catch.hpp"
#include "interaction.hpp"
#include "math_utility.hpp"



using namespace _pyves;



TEST_CASE("Interaction Test")
{
    {
        Particle p1(CARTESIAN::Zero(), CARTESIAN(-1,1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        Particle p2(CARTESIAN(1.5*nth_root<6>(2.f), 0, 0), CARTESIAN(1,1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        p1.gamma = PI/180*45;
        p2.gamma = PI/180*45;
        p1.kappa = 1;
        p2.kappa = 1;
        p1.epsilon = 1;
        p2.epsilon = 1;
        p1.sigma = 1;
        p2.sigma = 2;
        auto box = Box<PBC::ON>(10,10,10);
        auto e = interaction::potentialEnergy(p1, p2, box, 3);
        CHECK(e == Approx(-1));
    }
    {
        Particle p1(CARTESIAN::Zero(), CARTESIAN(-1,1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        Particle p2(CARTESIAN(1.5*nth_root<6>(2.f)+10, 0, 0), CARTESIAN(1,1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        p1.gamma = PI/180*45;
        p2.gamma = PI/180*45;
        p1.kappa = 1;
        p2.kappa = 1;
        p1.epsilon = 1;
        p2.epsilon = 1;
        p1.sigma = 1;
        p2.sigma = 2;
        auto box = Box<PBC::ON>(10,10,10);
        auto e = interaction::potentialEnergy(p1, p2, box, 3);
        CHECK(e == Approx(-1));
    }
    {
        Particle p1(CARTESIAN::Zero(), CARTESIAN(-1,-1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        Particle p2(CARTESIAN(1.5*nth_root<6>(2.f), 0, 0), CARTESIAN(1,-1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        p1.gamma = PI/180*45;
        p2.gamma = PI/180*45;
        p1.kappa = 1;
        p2.kappa = 1;
        p1.epsilon = 1;
        p2.epsilon = 1;
        p1.sigma = 1;
        p2.sigma = 2;
        auto box = Box<PBC::ON>(10,10,10);
        auto e = interaction::potentialEnergy(p1, p2, box, 3);
        CHECK(e == Approx(-1));
    }
    {
        Particle p1(CARTESIAN::Zero(), CARTESIAN(-1,-1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        Particle p2(CARTESIAN(1.5*nth_root<6>(2.f)+10, 0, 0), CARTESIAN(1,-1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        p1.gamma = PI/180*45;
        p2.gamma = PI/180*45;
        p1.kappa = 1;
        p2.kappa = 1;
        p1.epsilon = 1;
        p2.epsilon = 1;
        p1.sigma = 1;
        p2.sigma = 2;
        auto box = Box<PBC::ON>(10,10,10);
        auto e = interaction::potentialEnergy(p1, p2, box, 3);
        CHECK(e == Approx(-1));
    }
    {
        Particle p1(CARTESIAN::Zero(), CARTESIAN(0,1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        Particle p2(CARTESIAN(1*nth_root<6>(2.f), 0, 0), CARTESIAN(1,1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        p1.gamma = PI/180*0;
        p2.gamma = PI/180*45;
        p1.kappa = 1;
        p2.kappa = 1;
        p1.epsilon = 1;
        p2.epsilon = 1;
        p1.sigma = 1;
        p2.sigma = 1;
        auto box = Box<PBC::ON>(10,10,10);
        auto e = interaction::potentialEnergy(p1, p2, box, 3);
        CHECK(e == Approx(-1));
    }
    {
        Particle p1(CARTESIAN::Zero(), CARTESIAN(0,-1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        Particle p2(CARTESIAN(1*nth_root<6>(2.f), 0, 0), CARTESIAN(1,-1,0).normalized(), 1, 1, 1, 1, "UNDEF");
        p1.gamma = PI/180*0;
        p2.gamma = PI/180*45;
        p1.kappa = 1;
        p2.kappa = 1;
        p1.epsilon = 1;
        p2.epsilon = 1;
        p1.sigma = 1;
        p2.sigma = 1;
        auto box = Box<PBC::ON>(10,10,10);
        auto e = interaction::potentialEnergy(p1, p2, box, 3);
        CHECK(e == Approx(-1));
    }
}
