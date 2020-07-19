#define BOOST_TEST_MODULE Monomial Test
#include "boost/test/included/unit_test.hpp"
#include "boost/test/tools/output_test_stream.hpp"

#include "piranha.h"
// this has to also include the functionality inplemented by the BaseTerm class
// this class is derived from BaseTerm

//nearly all tests are testing BaseTerm except:
// multiply_test
//

//	template <class Cf, class Key, char Separator, class Allocator>
//class Monomial : public BaseTerm<Cf, Key, Separator, Allocator, Monomial<Cf, Key, Separator, Allocator> >
using namespace piranha;
using boost::test_tools::output_test_stream;


// using ExpoVector as key otherwise we have to create intermediate classes
// This could represent a polynomial monomial with integer exponents and double coefficients
namespace {
    using MonoType = Monomial<double_cf, ExpoVector<int, 0>, '|', std::allocator<char>>;
}

/////////////////////////////////////////////////
// Monomial specific tests
/////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(defaultconstruction_test)
{
//#define PIRANHA_TERM_CTORS(TermName) \
//	explicit TermName(): ancestor() {} \

//	template <class ArgsTuple> \
//	explicit TermName(const std::string &str, const ArgsTuple &argsTuple): \
//			ancestor(str, argsTuple) {} \

//	explicit TermName(const CfType &c, const KeyType &t): ancestor(c, t) {} \

//	template <class Cf2, class ArgsTuple> \
//	explicit TermName(const TermName<Cf2, KeyType, Separator, Allocator> &term, const ArgsTuple &a): \
//			ancestor(term, a) {} \

//	template <class Cf2, class Key2> \
//	explicit TermName(const TermName<Cf2, Key2, Separator, Allocator> &term): \
//			ancestor(term) {}

   MonoType mono;
   BOOST_TEST(mono.cf == 0);
   BOOST_TEST(mono.key.size() == 0); // cf and key are public in BaseTerm. Is that a good idea.
}

BOOST_AUTO_TEST_CASE(stringconstruction_test)
{
    // these are the same tests as for VectorKey
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    BOOST_CHECK_THROW(MonoType mono(" -1.0| 2, -3", polyArgs),boost::bad_lexical_cast); // has invalid separator for key
    BOOST_REQUIRE_THROW(MonoType mon0(" -1.0 | 2| 3", polyArgs), value_error); // has seaprator for cf|key twice; invalid
    BOOST_CHECK_NO_THROW(MonoType mon1(" -1.0| 2; -3", polyArgs)); // correct

    MonoType mono("-1.0|2;-3", polyArgs);
    BOOST_TEST(mono.cf.get_value() == -1.0);
    BOOST_TEST(mono.key.size() == 2);
    BOOST_TEST(mono.key[0] == 2);
    BOOST_TEST(mono.key[1] == -3);
}

BOOST_AUTO_TEST_CASE(cfandkeyconstruction_test)
{
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);
    ExpoVector<int, 0> key("1;2", polyArgs);
    double_cf value00(-1.0, polySym);

    MonoType mono(value00, key);
    BOOST_TEST(mono.cf == -1.0);
    BOOST_TEST(mono.key.size() == 2);
    BOOST_TEST(mono.key[0] == 1);
    BOOST_TEST(mono.key[1] == 2);
}

BOOST_AUTO_TEST_CASE(copyconstruction_test)
{
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);
    
    MonoType mono00("-1.0|1;2", polyArgs);

    MonoType mono01(mono00);
    BOOST_TEST((mono01 == mono00)); //this also tests operator==()
}

BOOST_AUTO_TEST_CASE(copyconstructionargs_test)
{
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    MonoType mono00("-1.0|1;2", polyArgs);

    MonoType mono01(mono00, polyArgs); // what is that actually good for???
    BOOST_TEST((mono01 == mono00)); //this also tests operator==()
}

BOOST_AUTO_TEST_CASE(mulitply_test)
{
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    MonoType mono00("-1.0|1;2", polyArgs);
    MonoType mono01("2.0|-2;1", polyArgs);

    MonoType::multiplication_result resultTuple; // this is atuple. here it is 1 1-tuple. why so complicated?? Anything else planned??

    // why not an operator ?????
    MonoType::multiply(mono00, mono01, resultTuple, polyArgs);

    MonoType result = resultTuple.get<0>();
    BOOST_TEST(result.cf == -2.0);
    BOOST_TEST(result.key.size() == 2);
    BOOST_TEST(result.key[0] == -1);
    BOOST_TEST(result.key[1] == 3);
}

////////////////////////////////////////////////////////////////////
// General tests. Covers basically BaseTerm
////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(get_test)
{
    // this get looks the same as if it would be for a tuple. Is there a purpose for this??
    // are actually both types of get needed????
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    MonoType mono00("-1.0|1;2", polyArgs);

    double_cf value = mono00.get<0>(); // get coefficient via echelon level. only 0 and 1 are supported. What is that good for????
    BOOST_TEST(value.get_value() == -1);
}

BOOST_AUTO_TEST_CASE(constget_test)
{
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    MonoType mono00("-1.0|1;2", polyArgs);

    double_cf const value = mono00.get<0>(); // get coefficient via echelon level. only 0 and 1 are supported. What is that good for????
    BOOST_TEST(value.get_value() == -1);
}

BOOST_AUTO_TEST_CASE(swap_test)
{
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    MonoType mono00("-1.0|1;2", polyArgs);
    MonoType mono01("2.0|2;3", polyArgs);
    MonoType const savedmono00(mono00);
    MonoType const savedmono01(mono01);

    mono00.swap(mono01);
    BOOST_TEST((mono01 == savedmono00)); // equality only covers the key!! why???
    BOOST_TEST(mono01.cf == savedmono00.cf);
    BOOST_TEST((mono00 == savedmono01));
    BOOST_TEST(mono00.cf == savedmono01.cf);
}

BOOST_AUTO_TEST_CASE(printPlain_test)
{
    Psym const e1("e1");
    Psym const e2("e2");
    Psym const e3("e3");
    VectorPsym const polySym = { e1, e2, e3 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    MonoType mono00("-1.0|2;3;4", polyArgs);
    output_test_stream output;
    BOOST_CHECK_NO_THROW(mono00.printPlain(output, polyArgs));
    BOOST_TEST(output.is_equal("-1|2;3;4"));
}

BOOST_AUTO_TEST_CASE(printPretty_test)
{
    Psym const e1("e1");
    Psym const e2("e2");
    Psym const e3("e3");
    VectorPsym const polySym = { e1, e2, e3 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    // unity for exponents
    MonoType mono00("-1.0|0;0;0", polyArgs);
    output_test_stream output;
    BOOST_CHECK_NO_THROW(mono00.printPretty(output, polyArgs));
    BOOST_TEST(output.is_equal("-1", true)); //empty output

    // drop cf because it is "1"
    MonoType mono01("1.0|2;3;4", polyArgs);
    BOOST_CHECK_NO_THROW(mono01.printPretty(output, polyArgs));
    BOOST_TEST(output.is_equal("e1**2*e2**3*e3**4", true));

    // prefix with "-" because of "-1.0"
    MonoType mono02("-1.0|2;3;4", polyArgs);
    BOOST_CHECK_NO_THROW(mono02.printPretty(output, polyArgs));
    BOOST_TEST(output.is_equal("-e1**2*e2**3*e3**4", true));

    // print everything
    MonoType mono03("-2.0|2;3;4", polyArgs);
    BOOST_CHECK_NO_THROW(mono03.printPretty(output, polyArgs));
    BOOST_TEST(output.is_equal("-2*e1**2*e2**3*e3**4", true));
}

BOOST_AUTO_TEST_CASE(printTex_test)
{
    Psym const e1("e1");
    Psym const e2("e2");
    Psym const e3("e3");
    VectorPsym const polySym = { e1, e2, e3 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    // unity for exponents
    MonoType mono00("-1.0|0;0;0", polyArgs);
    output_test_stream output;
    BOOST_CHECK_NO_THROW(mono00.printTex(output, polyArgs));
    BOOST_TEST(output.is_equal(" -1 ", true)); //just a blank prefixed! empty the output

    // drop cf because it is "1"
    MonoType mono01("1.0|2;3;4", polyArgs);
    BOOST_CHECK_NO_THROW(mono01.printTex(output, polyArgs));
    BOOST_TEST(output.is_equal(" e1 ^{2} e2 ^{3} e3 ^{4}", true)); // are the blanks really necessary?

    // prefix with "-" because of "-1.0"
    MonoType mono02("-1.0|2;3;4", polyArgs);
    BOOST_CHECK_NO_THROW(mono02.printTex(output, polyArgs));
    BOOST_TEST(output.is_equal("- e1 ^{2} e2 ^{3} e3 ^{4}", true));

    // print everything
    MonoType mono03("-2.0|2;3;4", polyArgs);
    BOOST_CHECK_NO_THROW(mono03.printTex(output, polyArgs));
    BOOST_TEST(output.is_equal(" -2  e1 ^{2} e2 ^{3} e3 ^{4}", true));
}

BOOST_AUTO_TEST_CASE(equality_test)
{
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    MonoType mono00("-1.0|1;2", polyArgs);
    MonoType mono01("2.0|1;2", polyArgs);

    BOOST_TEST((mono00 == mono01)); // equality operates on the keys on;y

}

BOOST_AUTO_TEST_CASE(isCanonical_test)
{
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    MonoType mono00("-1.0|1;2", polyArgs);
    BOOST_TEST(mono00.isCanonical(polyArgs)); // is always tru for BaseTerm == Monomial

}

BOOST_AUTO_TEST_CASE(iscanonicalise_test)
{

}

BOOST_AUTO_TEST_CASE(hasher_test)
{
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    MonoType mono00("-1.0|1;2", polyArgs);

    MonoType::Hasher hasher;
    BOOST_TEST((hasher(mono00) == hash_value(mono00)));
}

BOOST_AUTO_TEST_CASE(hash_test)
{
    Psym const e1("e1");
    Psym const e2("e2");

    VectorPsym const polySym = { e1, e2 };
    boost::tuple<VectorPsym> const polyArgs(polySym);

    MonoType mono00("-1.0|1;2", polyArgs);

    std::size_t const hash = hash_value(mono00); // namespace overload of metthod. What is this good for
}
