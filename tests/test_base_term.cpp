#define BOOST_TEST_MODULE base_term Test
#include "boost/test/included/unit_test.hpp"

// BaseTerm can not be directly tested becaue of CRTP
// see tests for derived versions of it that shall also cover 
// the tests in the base template
// Derived  are: FourierSeriesTerm
//               Monomial

// Should we create a suit eand put all the derived classes in there?
BOOST_AUTO_TEST_CASE(construction_test)
{
	BOOST_TEST_MESSAGE("See tests for derived classes");
}