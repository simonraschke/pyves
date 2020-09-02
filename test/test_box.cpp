#include "catch.hpp"
#include "box.hpp"



using namespace _pyves;



TEST_CASE("Box Constructor Test")
{
    Box<PBC::ON> box1;
    CHECK(box1.getLengthX()== Approx(static_cast<REAL>(0)));
    CHECK(box1.getLengthY()== Approx(static_cast<REAL>(0)));
    CHECK(box1.getLengthZ() == Approx(static_cast<REAL>(0)));

    Box<PBC::ON> box2(20,10,15);
    CHECK(box2.getLengthX()== Approx(static_cast<REAL>(20)));
    CHECK(box2.getLengthY()== Approx(static_cast<REAL>(10)));
    CHECK(box2.getLengthZ() == Approx(static_cast<REAL>(15)));

    Box<PBC::ON> box3(CARTESIAN(10,15,20));
    CHECK(box3.getLengthX()== Approx(static_cast<REAL>(10)));
    CHECK(box3.getLengthY()== Approx(static_cast<REAL>(15)));
    CHECK(box3.getLengthZ() == Approx(static_cast<REAL>(20)));
}



TEST_CASE("Box Simple Method Test")
{
    Box<PBC::ON> box1;
    Box<PBC::ON> box2(20,10,15);
    Box<PBC::ON> box3(CARTESIAN(10,15,20));

    box1.setLengthX(7);
    CHECK(box1.getLengthX() == Approx(static_cast<REAL>(7)));

    box1.setLengthY(3);
    CHECK(box1.getLengthY() == Approx(static_cast<REAL>(3)));

    box1.setLengthZ(6);
    CHECK(box1.getLengthZ() == Approx(static_cast<REAL>(6)));

    CHECK(box1.getVolume() == Approx(static_cast<REAL>(7*3*6)));
    CHECK(box2.getVolume() == Approx(static_cast<REAL>(20*10*15)));
    CHECK(box3.getVolume() == Approx(static_cast<REAL>(10*15*20)));
    
    CHECK((box2.getCenter()-CARTESIAN(10,5,7.5)).cwiseAbs().isApprox(CARTESIAN::Zero()));
}



TEST_CASE("Box PBC / NoPBC Test")
{
    Box<PBC::ON> boxPBC(10,10,10);
    Box<PBC::OFF> boxNoPBC(CARTESIAN(10,10,10));

    CHECK(boxPBC.distanceVector(CARTESIAN(0,0,0), CARTESIAN(1,1,1)).cwiseAbs().isApprox(CARTESIAN(1,1,1)));
    CHECK(boxPBC.distanceVector(CARTESIAN(0,0,0), CARTESIAN(11,11,11)).cwiseAbs().isApprox(CARTESIAN(1,1,1)));
    CHECK(boxNoPBC.distanceVector(CARTESIAN(0,0,0), CARTESIAN(1,1,1)).cwiseAbs().isApprox(CARTESIAN(1,1,1)));
    CHECK(boxNoPBC.distanceVector(CARTESIAN(0,0,0), CARTESIAN(11,11,11)).cwiseAbs().isApprox(CARTESIAN(11,11,11)));

    CHECK(boxPBC.distance(CARTESIAN(0,0,0), CARTESIAN(1,1,1)) == Approx(1.7320508075688772));
    CHECK(boxPBC.distance(CARTESIAN(0,0,0), CARTESIAN(11,11,11)) == Approx(1.7320508075688772));
    CHECK(boxNoPBC.distance(CARTESIAN(0,0,0), CARTESIAN(1,1,1)) == Approx(1.7320508075688772));
    CHECK(boxNoPBC.distance(CARTESIAN(0,0,0), CARTESIAN(11,11,11)) == Approx(19.05255888325765));

    CHECK(boxPBC.squaredDistance(CARTESIAN(0,0,0), CARTESIAN(1,1,1)) == Approx(3));
    CHECK(boxPBC.squaredDistance(CARTESIAN(0,0,0), CARTESIAN(11,11,11)) == Approx(3));
    CHECK(boxNoPBC.squaredDistance(CARTESIAN(0,0,0), CARTESIAN(1,1,1)) == Approx(3));
    CHECK(boxNoPBC.squaredDistance(CARTESIAN(0,0,0), CARTESIAN(11,11,11)) == Approx(363));

    CHECK(boxPBC.scaleToBox(CARTESIAN(11,11,11)).cwiseAbs().isApprox(CARTESIAN(1,1,1)));

    CHECK(boxPBC.contains(CARTESIAN(4,4,4)));
    CHECK(boxPBC.contains(CARTESIAN(10,10,10)));
    CHECK(boxPBC.contains(CARTESIAN(11,11,11)));
    CHECK(boxPBC.contains(CARTESIAN(1,-1,-21)));
    CHECK(boxNoPBC.contains(CARTESIAN(4,4,4)));
    CHECK(boxNoPBC.contains(CARTESIAN(10,10,10)));
    CHECK_FALSE(boxNoPBC.contains(CARTESIAN(11,11,11)));
    CHECK_FALSE(boxNoPBC.contains(CARTESIAN(1,-1,-21)));
    
    CHECK(boxNoPBC.contains(boxNoPBC.randomPointInside()));
}