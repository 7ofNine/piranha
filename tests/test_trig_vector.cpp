#define BOOST_TEST_MODULE trig_vector Test
#include "boost/test/included/unit_test.hpp"

#include "piranha.h"

using namespace std;
using namespace piranha;

namespace {

	using BaseTrigVector = TrigVector<int, 1>;
}


BOOST_AUTO_TEST_CASE(construction_test)
{
	BaseTrigVector constructed;
    BOOST_TEST(constructed.getFlavour() == true);
    BOOST_TEST(constructed.size() == 0);
}

BOOST_AUTO_TEST_CASE(construction_different_position_test) // is that actually used anywhater?
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(construct_from_string_list_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(construct_from_string_list2_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(construct_from_psym_test)
{
    // where is that very restricted construction used.
    // there is no explicit connection between the psym and the value 1 (its position!) on the container, only via external coding. dangerous
    Psym e1{"e1"};  //a exponential symbol
    Psym t1("t1");  // trigonometric symbol
    Psym t2("t2");  // trigonometric symbol

    VectorPsym trigSym = { t1 };
    VectorPsym expoSym = { e1 };

    //NTuple<VectorPsym, 2> expoTrigArgsType; // how to initialize a NTuple??
    //expoTrigArgsType expoTrigArgs = expoTrigargs.makeTuple(expoSym, trigSym);
    // expo keys are level 0, trig args are level 1
    boost::tuple<VectorPsym, VectorPsym> expoTrigArgs{ expoSym, trigSym };

    {
        // Is that really used anywhere ?
        BOOST_CHECK_NO_THROW( BaseTrigVector test(t1, 1, expoTrigArgs));

        // wrong position (echelon level) creates an empty key. Should we even allow this???
        BaseTrigVector tmp(t1, 0, expoTrigArgs);
        BOOST_TEST(tmp.size() == 0);

        // trig argument in echelon level 1
        BaseTrigVector tmp1(t1, 1, expoTrigArgs);
        BOOST_TEST(tmp1.size() == 1);
        BOOST_TEST(tmp1[0] == 1);
        BOOST_TEST(tmp1.getFlavour() == true);
    }

    // symbol doesn't match
    BOOST_REQUIRE_THROW(BaseTrigVector tmp(t2, 1, expoTrigArgs), assertion_error);

    {
        //position/echelon level doesn't match, acts like an empty constructor
        BaseTrigVector tmp(t1, 0, expoTrigArgs);
        BOOST_TEST(tmp.size() == 0);
    }

    // more than one symbol
    {
        VectorPsym trigSyms{ t1, t2 };
        boost::tuple<VectorPsym, VectorPsym> expoTrigArgs{ expoSym, trigSym };
        BOOST_REQUIRE_THROW(BaseTrigVector tmp(t2, 1, expoTrigArgs), assertion_error);
    }
}

BOOST_AUTO_TEST_CASE(multiply_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(flavour_test)
{
    BaseTrigVector constructed;
    // default is flavour (= cos)
    BOOST_TEST(constructed.getFlavour() == true);

    // we can reset/set it
    constructed.setFlavour(false);
    BOOST_TEST(constructed.getFlavour() == false);

    constructed.setFlavour(true);
    BOOST_TEST(constructed.getFlavour() == true);
}

BOOST_AUTO_TEST_CASE(printPlain_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(printPlainSorted_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(printTex_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(unity_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(harmonicDegree_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(partialHarmonicDegree_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(harmonicOrder_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(partialHarmonicOrder_test)
{
	BOOST_TEST_FAIL("No test implemented");
}


BOOST_AUTO_TEST_CASE(norm_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(signum_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(swap_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(trim_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(ignorable_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(comparison_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(hash_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(partialDerivative_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(power_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(substitution_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(eiSubstitution_test)
{
	BOOST_TEST_FAIL("No test implemented");
}
