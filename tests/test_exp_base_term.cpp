#define BOOST_TEST_MODULE BaseTerm Test
#include "boost/test/included/unit_test.hpp"

// Exp BaseTerm is not directly instantiated but is used as a base class for the derived classes
// using CRTP. I.e. the derived classes have to fulfill some interface requirements in order to be usable
// /derivable from a BaseTerm
// see tests for derived versions of it that shall also cover 
// the tests in the base template
// Derived  are: FourierSeriesTerm
//               Monomial
// 
// It makes more sense to test it's functions with the derived classes (?)

// this is for first experiments with using mp++ instead of mpir/mpfr directly



#include "piranha.h"
#include "mp++/mp++.hpp"
#include "exp_base_term.h"

using namespace std;
using namespace piranha;


namespace {
	class Mock {};

	using TestTerm = ExpBaseTerm< mppp::real, ExpoVector<int16_t, 0>, Mock>;
}
BOOST_AUTO_TEST_CASE(construction_test)
{
	
	TestTerm first;
}