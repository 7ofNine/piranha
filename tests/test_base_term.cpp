#define BOOST_TEST_MODULE BaseTerm Test
#include "boost/test/included/unit_test.hpp"

// BaseTerm is not directly instantiated but is used as a base class for the derived classes
// using CRTP. I.e. the derived classes have to fulfill some interface requirements in order to be usable
// /derivable from a BaseTerm
// see tests for derived versions of it that shall also cover 
// the tests in the base template
// Derived  are: FourierSeriesTerm
//               Monomial
// 
// It makes more sense to test it's functions with the derived classes (?)

#include "piranha.h"

using namespace std;
using namespace piranha;

BOOST_AUTO_TEST_CASE(construction_test)
{
	BOOST_TEST_MESSAGE("See tests for derived classes");
}