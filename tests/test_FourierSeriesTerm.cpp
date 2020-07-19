#define BOOST_TEST_MODULE FourierSeriesTerm Test
#include "boost/test/included/unit_test.hpp"
//#include "boost/test/tools/output_test_stream.hpp"


//
// Test for FourierSeriesTerms. BaseTerm is tested with Monomial
// Do we need addional tests for BaseSeries because we are using VectorKey insatead of ExpoKey??
// 

#include "piranha.h"

using namespace piranha;

namespace {
    using FourierType = FourierSeriesTerm<double_cf, TrigVector<int, 0>, '|', std::allocator<char>>; // we keep the echelon level at 0. i.e. just numerical coefficients
}

BOOST_AUTO_TEST_CASE(defaultconstruction_test)
{
    FourierType four;
    BOOST_TEST(four.cf == 0);
    BOOST_TEST(four.key.size() == 0); // cf and key are public in BaseTerm. Is that a good idea.
}

BOOST_AUTO_TEST_CASE(stringconstruction_test)
{
    // these are the same tests as for VectorKey
    Psym const t1("t1");
    Psym const t2("t2");

    VectorPsym const trigSym = { t1, t2 };
    boost::tuple<VectorPsym> const trigArgs(trigSym);

    BOOST_CHECK_THROW(FourierType four00(" -1.0| 2, -3", trigArgs), assertion_error); // has invalid separator for key and no flavour
    BOOST_CHECK_THROW(FourierType four01(" -1.0| 2, -3;s", trigArgs), assertion_error); // has invalid separator for key
    BOOST_CHECK_THROW(FourierType four02(" -1.0| 2; -3;d", trigArgs), value_error); // has invalid flavour type, actually tests TrigKey

    FourierType four03("-1.0|2;-3;s", trigArgs);
    BOOST_TEST(four03.cf.get_value() == -1.0);
    BOOST_TEST(four03.key.size() == 2);
    BOOST_TEST(four03.key[0] == 2);
    BOOST_TEST(four03.key[1] == -3);
    BOOST_TEST(four03.key.getFlavour() == false); // sin
}

BOOST_AUTO_TEST_CASE(cfandkeyconstruction_test)
{
    Psym const t1("t1");
    Psym const t2("t2");

    VectorPsym const trigSym = { t1, t2 };
    boost::tuple<VectorPsym> const trigArgs(trigSym);
    TrigVector<int, 0> key("1;2;s", trigArgs);
    double_cf value00(-1.0, trigSym);

    FourierType four(value00, key);
    BOOST_TEST(four.cf == -1.0);
    BOOST_TEST(four.key.size() == 2);
    BOOST_TEST(four.key[0] == 1);
    BOOST_TEST(four.key[1] == 2);
    BOOST_TEST(four.key.getFlavour() == false);
}


BOOST_AUTO_TEST_CASE(copyconstruction_test)
{
    Psym const t1("t1");
    Psym const t2("t2");

    VectorPsym const trigSym = { t1, t2 };
    boost::tuple<VectorPsym> const trigArgs(trigSym);

    FourierType four00("-1.0|1;2;s", trigArgs);

    FourierType four01(four00);
    BOOST_TEST((four01 == four00)); //this also tests operator==() but tests only the key
    BOOST_TEST((four00.cf == four01.cf));
}

BOOST_AUTO_TEST_CASE(copyconstructionarg_test)
{
    Psym const t1("t1");
    Psym const t2("t2");

    VectorPsym const trigSym = { t1, t2 };
    boost::tuple<VectorPsym> const trigArgs(trigSym);

    FourierType four00("-1.0|1;2;s", trigArgs);

    FourierType four01(four00, trigArgs);
    BOOST_TEST((four01 == four00)); //this also tests operator==() but tests only the key
    BOOST_TEST((four00.cf == four01.cf));
}

BOOST_AUTO_TEST_CASE(iscanonical_test)
{
    Psym const t1("t1");
    Psym const t2("t2");

    VectorPsym const trigSym = { t1, t2 };
    boost::tuple<VectorPsym> const trigArgs(trigSym);              

    FourierType four00("-1.0|1;2;s", trigArgs);                  // what is trigArgs good for???
    BOOST_TEST((four00.isCanonical(trigArgs) == true));

    FourierType four01("-1.0|-1;2;s", trigArgs);
    BOOST_TEST((four01.isCanonical(trigArgs) == false));

    FourierType four02("-1.0|1;-2;s", trigArgs);
    BOOST_TEST((four02.isCanonical(trigArgs) == true));
}

BOOST_AUTO_TEST_CASE(canonicalise_test)
{
    Psym const t1("t1");
    Psym const t2("t2");

    VectorPsym const trigSym = { t1, t2 };
    boost::tuple<VectorPsym> const trigArgs(trigSym);

    FourierType four00("-1.0|1;2;s", trigArgs);                  // what is trigArgs good for???
    four00.canonicalise(trigArgs);
    BOOST_TEST((four00.isCanonical(trigArgs) == true));          // nothing should have happened
    BOOST_TEST((four00.cf == -1.0));
    BOOST_TEST((four00.key[0] == 1));
    BOOST_TEST((four00.key[1] == 2));
    BOOST_TEST((four00.key.getFlavour() == false));

    FourierType four01("-1.0|-1;2;s", trigArgs);                  // what is trigArgs good for???
    four01.canonicalise(trigArgs);
    BOOST_TEST((four01.isCanonical(trigArgs) == true));
    BOOST_TEST((four01.cf == 1.0));                               // signum pulls through 
    BOOST_TEST(four01.key[0] == 1);
    BOOST_TEST(four01.key[1] == -2);
    BOOST_TEST((four01.key.getFlavour() == false)); // still sin

    FourierType four02("-1.0|-1;2;c", trigArgs);                  // what is trigArgs good for???
    four02.canonicalise(trigArgs);
    BOOST_TEST((four02.isCanonical(trigArgs) == true));
    BOOST_TEST((four02.cf == -1.0));                               // only the trigargs change their sign 
    BOOST_TEST(four02.key[0] == 1);
    BOOST_TEST(four02.key[1] == -2);
    BOOST_TEST((four02.key.getFlavour() == true)); // still sin

    FourierType four03("-1.0|1;2;c", trigArgs);                  // what is trigArgs good for???
    four03.canonicalise(trigArgs);
    BOOST_TEST((four03.isCanonical(trigArgs) == true));         //nothing has changed
    BOOST_TEST((four03.cf == -1.0));
    BOOST_TEST((four03.key[0] == 1));
    BOOST_TEST((four03.key[1] == 2));
    BOOST_TEST((four03.key.getFlavour() == true));
}

BOOST_AUTO_TEST_CASE(multiply_test)
{
    Psym const t1("t1");
    Psym const t2("t2");

    VectorPsym const trigSym = { t1, t2 };
    boost::tuple<VectorPsym> const trigArgs(trigSym);

    FourierType four00("-1|1;2;c", trigArgs);
    FourierType four01(" 1|2;3;s", trigArgs); // not the same in order to triger canonalization

    // the sequence of the results is not apriori defined it is implementation dependend
    // all two summands of the result are canonical and the one with the subtracted arguments comes first
    // this gives the result an order. The same way as I ormally remember kepping the subtraction at the second summand
    FourierType::multiplication_result result1;

    FourierType::multiply(four00, four00, result1, trigArgs);  // 1/2*(cos(a-b) + cos(a+b))
    FourierType result10 = result1.get<0>();
    FourierType result11 = result1.get<1>();
    BOOST_TEST(result10.isCanonical(trigArgs) == true);
    BOOST_TEST(result11.isCanonical(trigArgs) == true);
    BOOST_TEST(result10.cf == 0.5);
    BOOST_TEST(result11.cf == 0.5);
    BOOST_TEST(result10.key.getFlavour() == true);
    BOOST_TEST(result11.key.getFlavour() == true);
    BOOST_TEST(result10.key[0] == 0);    // a-b
    BOOST_TEST(result10.key[1] == 0);
    BOOST_TEST(result11.key[0] == 2);    // a+b
    BOOST_TEST(result11.key[1] == 4);

    FourierType::multiplication_result result2;
    FourierType::multiply(four00, four01, result2, trigArgs); // -1/2*sin(a-b) + 1/2*sin(a+b)
    FourierType result20 = result2.get<0>();
    FourierType result21 = result2.get<1>();
    BOOST_TEST(result20.isCanonical(trigArgs) == true);
    BOOST_TEST(result21.isCanonical(trigArgs) == true);
    BOOST_TEST(result20.cf == -0.5);
    BOOST_TEST(result21.cf == -0.5);
    BOOST_TEST(result20.key.getFlavour() == false);
    BOOST_TEST(result21.key.getFlavour() == false);
    BOOST_TEST(result20.key[0] == 1);    // a-b
    BOOST_TEST(result20.key[1] == 1);
    BOOST_TEST(result21.key[0] == 3);    // a+b
    BOOST_TEST(result21.key[1] == 5);

    
    FourierType::multiplication_result result3;
    FourierType::multiply(four01, four01, result3, trigArgs);  //1/2*cos(a-b) - 1/2*cos(a+b)
    FourierType result30 = result3.get<0>();
    FourierType result31 = result3.get<1>();
    BOOST_TEST(result30.isCanonical(trigArgs) == true);
    BOOST_TEST(result31.isCanonical(trigArgs) == true);
    BOOST_TEST(result30.cf == 0.5);
    BOOST_TEST(result31.cf == -0.5);
    BOOST_TEST(result30.key.getFlavour() == true);
    BOOST_TEST(result31.key.getFlavour() == true);
    BOOST_TEST(result30.key[0] == 0);    // a-b
    BOOST_TEST(result30.key[1] == 0);
    BOOST_TEST(result31.key[0] == 4);    // a+b
    BOOST_TEST(result31.key[1] == 6);

    FourierType::multiplication_result result4;   // this should be the same result as for result2 but not necessarily in the same sequence
    FourierType::multiply(four01, four00, result4, trigArgs); // 1/2*sin(a-b) + 1/2*sin(a+b)
    FourierType result40 = result4.get<0>();
    FourierType result41 = result4.get<1>();
    BOOST_TEST(result40.isCanonical(trigArgs) == true);
    BOOST_TEST(result41.isCanonical(trigArgs) == true);
    BOOST_TEST(result40.cf == -0.5);
    BOOST_TEST(result41.cf == -0.5);
    BOOST_TEST(result40.key.getFlavour() == false);
    BOOST_TEST(result41.key.getFlavour() == false);
    BOOST_TEST(result40.key[0] == 1);    // a-b
    BOOST_TEST(result40.key[1] == 1);
    BOOST_TEST(result41.key[0] == 3);    // a+b
    BOOST_TEST(result41.key[1] == 5);
}

