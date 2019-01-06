#define BOOST_TEST_MODULE BaseSeries Test
#include "boost/test/included/unit_test.hpp"

// This can not explicitly tested becasue of CRTP. The las parameter is the derived class parmeter
// This is probably incorrect. Create an intermediate class with just a bunch of definitons in it
BOOST_AUTO_TEST_CASE(construction_test)
{
	BOOST_TEST_FAIL("No test implemented");
}