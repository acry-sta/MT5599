#include <cxxtest/TestSuite.h>

class TestWritingTestsTutorial: public CxxTest::TestSuite
{
public:
    void TestOnePlusOneEqualsTwo()
    {
        int some_number = 1 + 1;
        TS_ASSERT_EQUALS(some_number, 2);
        double another_number = 1.000001 + 1.0001;
        TS_ASSERT_DELTA(another_number, 2.0, 1e-2);
    }
    void TestSomeOtherStuff()
    {
        TS_ASSERT(1==1); // however, it is better to use TS_ASSERT_EQUALS, below
        TS_ASSERT_EQUALS((true||false), true);
        TS_ASSERT_DIFFERS(1.348329534564385643543957436, 1.348329534564395643543957436);
        TS_ASSERT_LESS_THAN(2.71828183, 3.14159265); // Note: to test if x is greater than y, use TS_ASSERT_LESS_THAN(y,x)
        TS_ASSERT_LESS_THAN_EQUALS(-1e100, 1e100);
        TS_ASSERT_THROWS_ANYTHING(throw 0;); // normally you would put a function call inside the brackets

        unsigned x;
        // The following TS_ASSERT_THROWS_NOTHING may be useful if you want to be certain that there are no uncaught exceptions
        TS_ASSERT_THROWS_NOTHING(x=1;);  // normally you would put a function call inside the brackets
        TS_ASSERT_EQUALS(x, 1u); //Note that x and 1u are of the same type: unsigned integer
    }

    void donotTestThis()
    {
        TS_ASSERT_EQUALS(1,   2);
        TS_ASSERT_EQUALS(1u,  2u);
        TS_ASSERT_EQUALS(1.0, 2.0);
    }
};