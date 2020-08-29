#include "catch.hpp"
#include "point.hpp"



TEST_CASE("Point Constructor Test")
{
    SECTION("Normal Construction")
    {
        Point<int> p1(1,1,1);
        CHECK_NOTHROW(Point<int>());
        CHECK_NOTHROW(Point<int>(1,1,1));
        CHECK_NOTHROW(Point<int>(p1));

        Point<long> p2(1,1,1);
        CHECK_NOTHROW(Point<long>());
        CHECK_NOTHROW(Point<long>(1,1,1));
        CHECK_NOTHROW(Point<long>(p2));

        Point<float> p3(1.2, -1.3, 1.4);
        CHECK_NOTHROW(Point<float>());
        CHECK_NOTHROW(Point<float>(1.2, -1.3, 1.4));
        CHECK_NOTHROW(Point<float>(p3));

        Point<double> p4(1.2, -1.3, 1.4);
        CHECK_NOTHROW(Point<double>());
        CHECK_NOTHROW(Point<double>(1.2, -1.3, 1.4));
        CHECK_NOTHROW(Point<double>(p4));
    }

    SECTION("Conversions")
    {
        Point<int> p1(1,1,1);
        CHECK_NOTHROW(Point<int>(p1));
        CHECK_NOTHROW(Point<long>(p1));
        CHECK_NOTHROW(Point<float>(p1));
        CHECK_NOTHROW(Point<double>(p1));

        Point<long> p2(1,1,1);
        CHECK_NOTHROW(Point<int>(p2));
        CHECK_NOTHROW(Point<long>(p2));
        CHECK_NOTHROW(Point<float>(p2));
        CHECK_NOTHROW(Point<double>(p2));

        Point<float> p3(1.2, -1.3, 1.4);
        CHECK_NOTHROW(Point<int>(p3));
        CHECK_NOTHROW(Point<long>(p3));
        CHECK_NOTHROW(Point<float>(p3));
        CHECK_NOTHROW(Point<double>(p3));

        Point<double> p4(1.2, -1.3, 1.4);
        CHECK_NOTHROW(Point<int>(p4));
        CHECK_NOTHROW(Point<long>(p4));
        CHECK_NOTHROW(Point<float>(p4));
        CHECK_NOTHROW(Point<double>(p4));
    }
}



TEST_CASE("Point Operator Class Test")
{
    Point<int> p1(1,1,1);
    Point<int> p2(4,2,3);

    SECTION("operator+")
    {
        Point<int> p3 = p1 + p2;
        CHECK(p3.x == 5);
        CHECK(p3.y == 3);
        CHECK(p3.z == 4);
    }

    SECTION("operator-")
    {
        Point<int> p3 = p2 - p1;
        CHECK(p3.x == 3);
        CHECK(p3.y == 1);
        CHECK(p3.z == 2);
    }

    SECTION("operator+=")
    {
        Point<int> p3(0,2,1);
        Point<int> p4(3,7,-5);
        p3 += p4;
        CHECK(p3.x == 3);
        CHECK(p3.y == 9);
        CHECK(p3.z == -4);
    }
}



TEST_CASE("Point Member Function Test")
{
    Point<int> p1(4,2,3);
    CHECK(p1.dot(p1) == 29);
}