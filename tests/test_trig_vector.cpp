#define BOOST_TEST_MODULE trig_vector Test
#include "boost/test/included/unit_test.hpp"
#include "boost/test/output_test_stream.hpp"

#include "piranha.h"

using namespace std;
using namespace piranha;
using namespace piranha::manipulators;
using boost::test_tools::output_test_stream;

namespace {

	using BaseTrigVector = TrigVector<int, 1>;
    using HighTrigVector = TrigVector<int, 2>;
}


BOOST_AUTO_TEST_CASE(construction_test)
{
	BaseTrigVector constructed;
    BOOST_TEST(constructed.getFlavour() == true);
    BOOST_TEST(constructed.size() == 0);
}

BOOST_AUTO_TEST_CASE(construction_different_position_test) // is that actually used anywhater. 
{
    HighTrigVector temp2;
    temp2.resize(3);
    temp2[0] = 0;
    temp2[1] = 1;
    temp2[2] = 2;
    temp2.setFlavour(false);
    BaseTrigVector temp1(temp2);

    BOOST_TEST(temp1.size() == 3);
    BOOST_TEST(temp1[0] == 0);
    BOOST_TEST(temp1[1] == 1);
    BOOST_TEST(temp1[2] == 2);
    BOOST_TEST(temp1.getFlavour() == false);
}

BOOST_AUTO_TEST_CASE(construct_from_string_list_test)
{
    string const inlist("1; -6; c"); // separator in TrigVector is ";"> can we read it
    Psym const t1("t1");  // trigonometric symbol
    Psym const t2("t2");  // trigonometric symbol

    VectorPsym trigSym = { t1, t2 };

    boost::tuple<VectorPsym, VectorPsym> trigArgs{ trigSym, trigSym }; // if we don't have both levels the constructor compilation failes!
    BOOST_CHECK_NO_THROW(BaseTrigVector temp(inlist, trigArgs));
    BaseTrigVector temp(inlist, trigArgs);
    BOOST_CHECK(temp.size() == 2);
    BOOST_CHECK(temp[0] == 1);
    BOOST_CHECK(temp[1] == -6);
    BOOST_CHECK(temp.getFlavour() == true);
    
    string inlist2("1; -6; s");
    BaseTrigVector temp1(inlist2, trigArgs);
    BOOST_CHECK(temp1.size() == 2);
    BOOST_CHECK(temp1[0] == 1);
    BOOST_CHECK(temp1[1] == -6);
    BOOST_CHECK(temp1.getFlavour() == false);

    string inlist3("1; -6; x");
    BOOST_CHECK_THROW(BaseTrigVector temp(inlist3, trigArgs), value_error);

    string inlist4("1; -6; 7; c");
    BOOST_CHECK_THROW(BaseTrigVector temp(inlist4, trigArgs), assertion_error);

    // does that type of construction make sense anyhwere or is it even an error? investigate.
    // we need empty trigArgs, when would that happen? error by nature??
    //BaseTrigVector tempe("", trigArgs);
    //BOOST_CHECK(tempe.size() == 0);
    //BOOST_CHECK(temp.getFlavour() == true);

}

BOOST_AUTO_TEST_CASE(construct_from_string_list2_test)
{
    // this is nearly the same as construct_from_string_list_test but without argsTuple
    string inlist("1; -6; c"); // separator in TrigVector is ";"> can we read it, the blank before the flavour is a test!

    BaseTrigVector temp(inlist);

    BOOST_CHECK(temp.size() == 2);
    BOOST_CHECK(temp[0] == 1);
    BOOST_CHECK(temp[1] == -6);
    BOOST_CHECK(temp.getFlavour() == true);

    string inlist1("1; -6; s"); // separator in TrigVector is ";"> can we read it, the blank before the flavour is a test!

    BaseTrigVector temp1(inlist1);

    BOOST_CHECK(temp1.size() == 2);
    BOOST_CHECK(temp1[0] == 1);
    BOOST_CHECK(temp1[1] == -6);
    BOOST_CHECK(temp1.getFlavour() == false);

    string inlistx("1; -6; x"); // separator in TrigVector is ";"> can we read it, the blank before the flavour is a test!

    BOOST_CHECK_THROW(BaseTrigVector tempx(inlistx), value_error);

    // the "wrong" number of elements is not detectable! No test for that 

    //BaseTrigVector tempe("");
    //BOOST_CHECK(tempe.size() == 0);
    //BOOST_CHECK(temp.getFlavour() == true);
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
    BaseTrigVector multiplicand;
    BaseTrigVector multiplier;
    multiplier.resize(1);
    BaseTrigVector result1;
    BaseTrigVector result2;

    BOOST_REQUIRE_THROW(multiplicand.multiply(multiplier, result1, result2), assertion_error);

    BaseTrigVector multiplicand1;
    multiplicand1.resize(2);
    multiplicand1[0] = 1;
    multiplicand1[1] = 2;
    BaseTrigVector multiplier1;
    multiplier1.resize(2);
    multiplier1[0] = -1;
    multiplier1[1] = -2;
    BaseTrigVector result11;
    BaseTrigVector result12;
    BOOST_CHECK_NO_THROW(multiplicand1.multiply(multiplier1, result11, result12));
    BOOST_TEST(result11.size() == 2);
    BOOST_TEST(result11[0] == 2);
    BOOST_TEST(result11[1] == 4);

    BOOST_TEST(result12.size() == 2);
    BOOST_TEST(result12[0] == 0);
    BOOST_TEST(result12[1] == 0);
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
    // corresponds to printElements in vectorKey
    BaseTrigVector temp;
    temp.resize(3);
    temp[0] = 1;
    temp[1] = 2;
    temp[2] = 3;

    // we have to have symbols
    Psym e1{"e1"};
    Psym t1("t1");
    Psym t2("t2");
    Psym t3("t3");
    VectorPsym polySym{ e1 };
    VectorPsym trigSym{ t1,t2,t3 };

    boost::tuple<VectorPsym, VectorPsym> polyTrigArgs(polySym, trigSym);

    output_test_stream output;
    BOOST_CHECK_NO_THROW(temp.printPlain(output, polyTrigArgs));
    BOOST_TEST(output.is_equal("1;2;3;c"));

    temp.invertSign();
    BOOST_CHECK_NO_THROW(temp.printPlain(output, polyTrigArgs));
    BOOST_TEST(output.is_equal("-1;-2;-3;c"));

    // wrong number of args, they have to agree with the number of values in the vector
    VectorPsym trigSym2{ t1,t3 };

    boost::tuple<VectorPsym, VectorPsym> polyTrig2Args(polySym, trigSym2);
    BOOST_CHECK_THROW(temp.printPlain(output, polyTrig2Args), assertion_error);

    VectorPsym trigSym4{ t1,t3,t3,t3 };
    boost::tuple<VectorPsym, VectorPsym> polyTrig4Args(polySym, trigSym4);
    BOOST_CHECK_THROW(temp.printPlain(output, polyTrig4Args), assertion_error);

}

BOOST_AUTO_TEST_CASE(printPlainSorted_test)
{
    BaseTrigVector temp;
    temp.resize(3);
    temp[0] = 1;
    temp[1] = 2;
    temp[2] = 3;
    Psym t1("t1");
    Psym t2("t2");
    Psym t3("t3");
    Psym e1("e1");
    VectorPsym trigSym{ t1, t2, t3 };
    VectorPsym polySym{ e1 };
    boost::tuple<VectorPsym, VectorPsym> polyTrigArgs(polySym, trigSym);

    output_test_stream output;
    std::vector<std::pair<bool, std::size_t>> sort = { {true, 2}, {true, 0}, {true, 1} };

    BOOST_CHECK_NO_THROW(temp.printPlainSorted(output, sort, polyTrigArgs));
    BOOST_TEST(!output.is_empty(false));
    BOOST_TEST(output.is_equal("     3      1      2 c"));


    //sizes don't agree
    BOOST_CHECK_THROW(temp.printPlainSorted(output, std::vector<std::pair<bool, std::size_t>>(), polyTrigArgs), assertion_error);

}

BOOST_AUTO_TEST_CASE(printPrettyx) // boost test doesn't like printPretty?
{
        Psym t1("t1");
        Psym t2("t2");
        Psym t3("t3");
        Psym e1("e1");
        VectorPsym trigSym = { t1, t3, t2 };
        VectorPsym manySym = { t1,t1,t1,t1 };
        VectorPsym noSym;
        VectorPsym polySym = { e1 };
        //boost::tuple<VectorPsym,VectorPsym> polyArgs(noSym, trigSym);
        boost::tuple<VectorPsym,VectorPsym> manyArgs(noSym, manySym);
        boost::tuple<VectorPsym,VectorPsym> trigArgs(polySym,trigSym);
        boost::tuple<VectorPsym,VectorPsym> noArgs(noSym, noSym);

        BaseTrigVector term;
        term.resize(3);
        term[0] = 1;
        term[1] = 1;
        term[2] = 1;
        output_test_stream output;

        // throw if sizes between ExpoVector and ArgsTuple component
        // are different
        BOOST_CHECK_THROW(term.printPretty(output, noArgs), assertion_error);
        BOOST_CHECK_THROW(term.printPretty(output, manyArgs), assertion_error);

        // everything is exponent 1
        BOOST_CHECK_NO_THROW(term.printPretty(output, trigArgs));
        BOOST_TEST(!output.is_empty(false));
        BOOST_TEST(output.is_equal("cos(t1+t3+t2)"));

        //// leading factor is 0
        term[0] = 0;
        term[1] = 1;
        term[2] = 1;
        BOOST_CHECK_NO_THROW(term.printPretty(output, trigArgs));
        BOOST_TEST(!output.is_empty(false));
        BOOST_TEST(output.is_equal("cos(t3+t2)"));

        // leading factor is 0 and last is >1
        term[0] = 0;
        term[1] = 1;
        term[2] = 3;
        BOOST_CHECK_NO_THROW(term.printPretty(output, trigArgs));
        BOOST_TEST(!output.is_empty(false));
        BOOST_TEST(output.is_equal("cos(t3+3*t2)"));

        // middle one is >1
        term[0] = 1;
        term[1] = 3;
        term[2] = 1;
        BOOST_CHECK_NO_THROW(term.printPretty(output, trigArgs));
        BOOST_TEST(!output.is_empty(false));
        BOOST_TEST(output.is_equal("cos(t1+3*t3+t2)"));

        // test a sin
        term.setFlavour(false);
        term[0] = 1;
        term[1] = 3;
        term[2] = 1;
        BOOST_CHECK_NO_THROW(term.printPretty(output, trigArgs));
        BOOST_TEST(!output.is_empty(false));
        BOOST_TEST(output.is_equal("sin(t1+3*t3+t2)"));

        // test a sin and negative factors
        term.setFlavour(false);
        term[0] = -1;
        term[1] = -3;
        term[2] = -1;
        BOOST_CHECK_NO_THROW(term.printPretty(output, trigArgs));
        BOOST_TEST(!output.is_empty(false));
        BOOST_TEST(output.is_equal("sin(-t1-3*t3-t2)"));

}


BOOST_AUTO_TEST_CASE(printTex_test)
{
    Psym t1("t1");
    Psym t2("t2");
    Psym t3("t3");
    VectorPsym trigSym = { t1, t2, t3 };
    VectorPsym manySym = { t1,t1,t1,t1 };
    VectorPsym noSym;
    boost::tuple<VectorPsym, VectorPsym> trigArgs(noSym, trigSym);
    boost::tuple<VectorPsym, VectorPsym> manyArgs(noSym, manySym);
    boost::tuple<VectorPsym, VectorPsym> noArgs(noSym, noSym);

    BaseTrigVector term;
    term.resize(3);
    term[0] = 1;
    term[1] = 1;
    term[2] = 1;
    output_test_stream output;

    // throw if sizes between TrigVector and ArgsTuple component
    // are different
    BOOST_CHECK_THROW(term.printTex(output, noArgs), assertion_error);
    BOOST_CHECK_THROW(term.printTex(output, manyArgs), assertion_error);

    // everything is exponent 1
    BOOST_CHECK_NO_THROW(term.printTex(output, trigArgs));
    BOOST_TEST(!output.is_empty(false));
    BOOST_TEST(output.is_equal("\\cos\\left(t1+t2+t3\\right)")); // escape sequences!

    //the same for sin;
    term.setFlavour(false);
    BOOST_CHECK_NO_THROW(term.printTex(output, trigArgs));
    BOOST_TEST(!output.is_empty(false));
    BOOST_TEST(output.is_equal("\\sin\\left(t1+t2+t3\\right)")); // escape sequences!

    //// leading factor is 0
    term[0] = 0;
    term[1] = 1;
    term[2] = 1;
    BOOST_CHECK_NO_THROW(term.printTex(output, trigArgs));
    BOOST_TEST(!output.is_empty(false));
    BOOST_TEST(output.is_equal("\\sin\\left(t2+t3\\right)"));

    // leading exponent is 0 and last is >1
    term[0] = 0;
    term[1] = 1;
    term[2] = 3;
    term.setFlavour(true); // back to cos
    BOOST_CHECK_NO_THROW(term.printTex(output, trigArgs));
    BOOST_TEST(!output.is_empty(false));
    BOOST_TEST(output.is_equal("\\cos\\left(t2+3t3\\right)"));

    // middle one is >1
    term[0] = 1;
    term[1] = 3;
    term[2] = 1;
    BOOST_CHECK_NO_THROW(term.printTex(output, trigArgs));
    BOOST_TEST(!output.is_empty(false));
    BOOST_TEST(output.is_equal("\\cos\\left(t1+3t2+t3\\right)"));

    // middle one is >1 and negative numbers
    term[0] = -1;
    term[1] = -3;
    term[2] = -1;
    BOOST_CHECK_NO_THROW(term.printTex(output, trigArgs));
    BOOST_TEST(!output.is_empty(false));
    BOOST_TEST(output.is_equal("\\cos\\left(-t1-3t2-t3\\right)"));
}

BOOST_AUTO_TEST_CASE(unity_test)
{
    BaseTrigVector temp;
    temp.resize(3);
    temp[0] = 1;
    temp[1] = 1;
    temp[2] = 1;
    BOOST_TEST(temp.isUnity() == false);

    temp[0] = 0;
    temp[1] = 0;
    temp[2] = 0;

    BOOST_TEST(temp.isUnity() == true);

    temp.setFlavour(false);
    BOOST_TEST(temp.isUnity() == false);

}

BOOST_AUTO_TEST_CASE(harmonicDegree_test)
{
    BaseTrigVector temp;
    temp.resize(3);
    temp[0] = 1;
    temp[1] = 2;
    temp[2] = 3;
    BOOST_TEST(temp.harmonicDegree() == 6);

    BaseTrigVector temp1;
    temp1.resize(3);
    temp1[0] = -1;
    temp1[1] = -2;
    temp1[2] = 3;
    BOOST_TEST(temp1.harmonicDegree() == 0);

    BaseTrigVector temp2;
    temp2.resize(3);
    temp2[0] = -1;
    temp2[1] = -2;
    temp2[2] = -3;
    BOOST_TEST(temp2.harmonicDegree() == -6);
}

BOOST_AUTO_TEST_CASE(partialHarmonicDegree_test)
{
    Psym t1("t1");
    Psym t2("t2");
    Psym t3("t3");

    VectorPsym trigSym{ t1, t2, t3 };
    VectorPsym noSym;
    boost::tuple<VectorPsym, VectorPsym> polyTrigArgs(noSym, trigSym);

    VectorPsym symbolsUsed{ t1, t3 };
    auto const posTuple = psyms2pos(symbolsUsed, polyTrigArgs);

    BaseTrigVector temp;
    temp.resize(3);
    temp[0] = 1;
    temp[1] = 2;
    temp[2] = 3;

	BOOST_TEST(temp.partialHarmonicDegree(posTuple)==4);
}

BOOST_AUTO_TEST_CASE(harmonicOrder_test)
{
    // harmonic order is defined to be the same has harmonic degree TODO:: why does it exist??

    BaseTrigVector temp;
    temp.resize(3);
    temp[0] = 1;
    temp[1] = 2;
    temp[2] = 3;
    BOOST_TEST(temp.harmonicDegree() == temp.harmonicOrder());

    BaseTrigVector temp1;
    temp1.resize(3);
    temp1[0] = -1;
    temp1[1] = -2;
    temp1[2] = 3;
    BOOST_TEST(temp.harmonicDegree() == temp.harmonicOrder());

    BaseTrigVector temp2;
    temp2.resize(3);
    temp2[0] = -1;
    temp2[1] = -2;
    temp2[2] = -3;
    BOOST_TEST((temp.harmonicDegree() == temp.harmonicOrder()));

}

BOOST_AUTO_TEST_CASE(partialHarmonicOrder_test)
{
    // harmonic order is defined to be the same has harmonic degree TODO:: why does it exist?? Only for the use with a toolbox
    Psym t1("t1");
    Psym t2("t2");
    Psym t3("t3");

    VectorPsym trigSym{ t1, t2, t3 };
    VectorPsym noSym;
    boost::tuple<VectorPsym, VectorPsym> polyTrigArgs(noSym, trigSym);

    VectorPsym symbolsUsed{ t1, t3 };
    auto const posTuple = psyms2pos(symbolsUsed, polyTrigArgs);

    BaseTrigVector temp;
    temp.resize(3);
    temp[0] = 1;
    temp[1] = 2;
    temp[2] = 3;

    BOOST_TEST(temp.partialHarmonicDegree(posTuple) == temp.partialHarmonicOrder(posTuple));

}


BOOST_AUTO_TEST_CASE(norm_test)
{
    // where is that very restricted construction used.
// there is no explicit connection between the psym and the value 1 (its position!) on the container, only via external coding. dangerous
    Psym e1{ "e1" };  //a exponential symbol
    Psym t1("t1");  // trigonometric symbol
    Psym t2("t2");  // trigonometric symbol

    VectorPsym trigSym = { t1, t2 };
    VectorPsym expoSym = { e1 };

    // expo keys are level 0, trig args are level 1
    boost::tuple<VectorPsym, VectorPsym> expoTrigArgs{ expoSym, trigSym };

    BaseTrigVector temp;
    temp.resize(3);
    BOOST_CHECK_THROW(temp.norm(expoTrigArgs), assertion_error);

    BaseTrigVector temp1;
    temp1.resize(1);
    temp[0] = 99;
    BOOST_TEST(temp1.norm(expoTrigArgs) == 1);
}

BOOST_AUTO_TEST_CASE(signum_test)
{
    BaseTrigVector test;
    test.resize(3);
    test[0] = 1;
    test[1] = -2;
    test[2] = 3;

    BaseTrigVector result(-test);
    BOOST_TEST(result.size() == test.size());
    BOOST_TEST(result[0] == -test[0]);
    BOOST_TEST(result[1] == -test[1]);
    BOOST_TEST(result[2] == -test[2]);
    BOOST_TEST(result.getFlavour() == test.getFlavour());

}

BOOST_AUTO_TEST_CASE(swap_test)
{
    BaseTrigVector vec1;
    vec1.resize(2);
    vec1[0] = 1;
    vec1[1] = -2;

    BaseTrigVector vec2;
    vec2.resize(2);
    vec2[0] = -1;
    vec2[1] = 2;
    vec2.setFlavour(false);
    
    BaseTrigVector const savedVec1(vec1);
    BaseTrigVector const savedVec2(vec2);

    vec1.swap(vec2);
    BOOST_TEST(vec1 == savedVec2);   // this also uses comparison
    BOOST_TEST(vec2 == savedVec1);
}

BOOST_AUTO_TEST_CASE(trim_test)
{
	BOOST_TEST_FAIL("No test implemented");
}

BOOST_AUTO_TEST_CASE(ignorable_test)
{
    Psym e1{ "e1" };  //a exponential symbol
    Psym t1("t1");  // trigonometric symbol
    Psym t2("t2");  // trigonometric symbol

    VectorPsym trigSym = { t1, t2 };
    VectorPsym expoSym = { e1 };

    // expo keys are level 0, trig args are level 1
    boost::tuple<VectorPsym, VectorPsym> expoTrigArgs{ expoSym, trigSym };

    BaseTrigVector ignore;
    ignore.resize(2);
    ignore[0] = 0;
    ignore[1] = 0;
    ignore.setFlavour(false);
    
    BOOST_TEST((ignore.isIgnorable(expoTrigArgs) == true)); // argsTuple is nowhere used in the method!!

    BaseTrigVector ignore1;
    ignore1.resize(2);
    ignore1[0] = 1;
    ignore1[1] = 0;
    ignore1.setFlavour(false);

    BOOST_TEST((ignore1.isIgnorable(expoTrigArgs) == false)); // argsTuple is nowhere used in the method!!

    BaseTrigVector ignore2;
    ignore2.resize(2);
    ignore2[0] = 0;
    ignore2[1] = 0;

    BOOST_TEST((ignore2.isIgnorable(expoTrigArgs) == false)); // argsTuple is nowhere used in the method!!

    BaseTrigVector ignore3;
    ignore3.resize(2);
    ignore3[0] = 1;
    ignore3[1] = 1;

    BOOST_TEST((ignore2.isIgnorable(expoTrigArgs) == false)); // argsTuple is nowhere used in the method!!
}

BOOST_AUTO_TEST_CASE(comparison_test)
{
    //equality
    BaseTrigVector vec1;
    vec1.resize(2);
    vec1[0] = 1;
    vec1[1] = -2;
    BaseTrigVector vecComp(vec1);
    BOOST_TEST(vec1 == vecComp);
    
    vecComp.setFlavour(false); // change to sin
    BOOST_TEST(((vec1 == vecComp) == false));

    BaseTrigVector vec2;
    vec2.resize(3);
    vec2[0] = 1;
    vec2[1] = -2;
    vec2[2] = 0;
    BOOST_TEST(((vec1 == vec2)==false));

    BaseTrigVector vec3;
    vec3.resize(2);
    vec3[0] = 1;
    vec3[1] = -3;
    BOOST_TEST(((vec1 == vec3) == false));

    // less
    BaseTrigVector vec4;
    vec4.resize(2);
    vec4[0] = 1;
    vec4[1] = -2;

    BaseTrigVector vec5;
    vec5.resize(2);
    vec5.setFlavour(false);
    vec5[0] = 1;
    vec5[1] = -2;

    BOOST_TEST(vec5 < vec4);
    BOOST_TEST(!(vec4 < vec5));

    BOOST_TEST(!(vec4 < vec4));
    BOOST_TEST(!(vec5 < vec5));

    BaseTrigVector vec6;
    vec6.resize(2);
    vec6.setFlavour(false);
    vec6[0] = -1;
    vec6[1] = -2;
    BOOST_TEST(vec6 < vec5);
}

BOOST_AUTO_TEST_CASE(hash_test)
{
    BaseTrigVector vec;
    vec.resize(3);
    vec[0] = -1;
    vec[1] = 2;
    vec[2] = -3;
    size_t value = vec.hash_value();  // just to execute the method
}









BOOST_AUTO_TEST_CASE(partialDerivative_test)
{
// TODO: In general is that a good idea that all this series handling is on the level of TrigVector???
// requires reconsdieraton but first cretae tests
//
// how to test this? it requires a series in the 
// call to partial as a template parameter.
// Is that in general a good construction
// that lower level classes require definitions from 
// higher levels?
// It requires Series as well as BaseSeries
// concepts?

// we use predefined qpoly. This also shows the problem with hierarchy
// PosTuple is used only from NAmedSeries it comes from method psym2Pos
// i.e. this method has many prerequisites.

    Psym e1("n1", "1;2;3");
    Psym e2("n2");
    Psym e3("n3", "4.0;5.0;6.0", 7);
    Psym e4("n4");
    Psym e5("n5");

    Psym t1("t1");
    Psym t2("t2");
    Psym t3("t3");
    Psym t4("t4");
    Psym t5("t5");

    Psym d1("d1");
    Psym d2("d2");
    Psym d3("d3");
    Psym d4("d4");
    Psym d5("d5");

    VectorPsym polySym{ e1, e3, e5 };
    VectorPsym trigSym{ t1, t3, t5 };
    VectorPsym divSym{ d1, d3, d5 };

    boost::tuple<VectorPsym> polyOnlyArgs(polySym);
    boost::tuple<VectorPsym, VectorPsym> polyTrigArgs(polySym, trigSym);
    boost::tuple<VectorPsym, VectorPsym, VectorPsym> polyTrigDivArgs(polySym, trigSym, divSym);
    boost::tuple<VectorPsym, VectorPsym, VectorPsym> allArgs{ {e1,e2,e3,e4,e5}, {t1,t2,t3,t4,t5}, {d1,d2,d3,d4,d5} };
    VectorPsym onePoly{ e3 };
    VectorPsym oneTrig{ t3 };

    auto  posTuple = psyms2pos(oneTrig, polyTrigDivArgs);

    qps result;
    qps::TermType::KeyType temp; // using BaseTrigVector leads to compilation error, use the TrigVector type declared in qps
                                // another indication of wrong dependency. This type is TrigVector<short, 1>
    temp.resize(3); // this is a cos term
    temp[0] = 3;
    temp[1] = 4;
    temp[2] = 5;

    result = temp.partial<qps>(posTuple, polyTrigDivArgs); // this result is actually a series. Is that the proper place for this test
                                                           // this is only partial differentiation after t3, why are we having the whole pos tuple in here
                                                           // instead of just the position we are interested in i.e. "t3" 

    BOOST_TEST(result.length() == 1); // should only have a single element length is from BaseSeries.
                                      // Most of the method should be from  BaseSeries.
    qps::TermType term;      // TermType is of class piranha::Monomial which is of type BaseTerm except for multiply method
    term = *(result.begin());  // get the only term
    qps::TermType::KeyType resultKey = term.get<1>(); // get the key <1> gives the key
    qps::TermType::CfType::TermType::CfType coefficient = (*(term.get<0>().begin())).get<0>();  // awkward
    BOOST_TEST((resultKey.getFlavour() == false)); 
    BOOST_TEST(coefficient == -4);
    BOOST_TEST(resultKey[0] == 3);
    BOOST_TEST(resultKey[1] == 4);
    BOOST_TEST(resultKey[2] == 5);

    qps::TermType::KeyType temp1; // using BaseTrigVector leads to compilation error, use the TrigVector type declared in qps
                                // another indication of wrong dependency. This type is TrigVector<short, 1>
    temp1.resize(3); // this is a cos term
    temp1[0] = 3;
    temp1[1] = 4;
    temp1[2] = 5;
    temp1.setFlavour(false);

    result = temp1.partial<qps>(posTuple, polyTrigDivArgs); // this result is actually a series. Is that the proper place for this test
    qps::TermType term1;      // TermType is of class piranha::Monomial which is of type BaseTerm except for multiply method
    term1 = *(result.begin());  // get the only term
    qps::TermType::KeyType resultKey1 = term1.get<1>(); // get the key <1> gives the key
    //qps::TermType::CfType::TermType::CfType test = (*(term1.get<0>().begin())).get<0>();
    qps::TermType::CfType::TermType::CfType coefficient1 = (*(term1.get<0>().begin())).get<0>();  // really akward. Only works for Poissonseries we should not return a series. That should be much higher
    BOOST_TEST((resultKey1.getFlavour() == true));
    BOOST_TEST(coefficient1 == 4); 
    BOOST_TEST(resultKey[0] == 3);
    BOOST_TEST(resultKey[1] == 4);
    BOOST_TEST(resultKey[2] == 5);


    VectorPsym twoTrig{ t3, t5 };
    auto  posTuple3 = psyms2pos(twoTrig, polyTrigDivArgs);
    qps result3;
    qps::TermType::KeyType temp3; // using BaseTrigVector leads to compilation error, use the TrigVector type declared in qps
                                // another indication of wrong dependency. This type is TrigVector<short, 1>
    temp3.resize(3); // this is a cos term
    temp3[0] = 3;
    temp3[1] = 4;
    temp3[2] = 5;

    BOOST_CHECK_THROW(result3 = temp3.partial<qps>(posTuple3, polyTrigDivArgs), assertion_error); // we only allow one variable to derive of

    VectorPsym lastTrig{ t5 };
    auto  posTuple4 = psyms2pos(lastTrig, polyTrigDivArgs);
    qps result4;
    qps::TermType::KeyType temp4; // using BaseTrigVector leads to compilation error, use the TrigVector type declared in qps
                                // another indication of wrong dependency. This type is TrigVector<short, 1>
    temp4.resize(1); // this is a cos term
    temp4[0] = 3;

    BOOST_CHECK_THROW(result4 = temp4.partial<qps>(posTuple4, polyTrigDivArgs), assertion_error); // outisde of range i.e. variable not in TrigVector
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

BOOST_AUTO_TEST_CASE(eval_test)
{
    BOOST_TEST_FAIL("No test implemented");
}