#define BOOST_TEST_MODULE Double coefficients
#include "boost/test/included/unit_test.hpp"

// contains also all test cases for the base class NumericalContainer
// NUmerical Container is only used as base class for other classes and can not be instantaiated by itself


#include <string>
#include "piranha.h"

using namespace piranha;
using namespace std;
//explicit double_cf() : ancestor() {}

//template <class T, class ArgsTuple>
//explicit double_cf(const T &x, const ArgsTuple &argsTuple) : ancestor(x, argsTuple) {}

//template <class ArgsTuple>
//explicit double_cf(const mp_rational &q, const ArgsTuple &argsTuple) : ancestor(q.to_double(), argsTuple) {}

//template <class ArgsTuple>
//explicit double_cf(const mp_integer &z, const ArgsTuple &argsTuple) : ancestor(z.to_double(), argsTuple) {}

//template <class ArgsTuple>
//explicit double_cf(const Psym &p, const int &n, const ArgsTuple &a) : ancestor(p, n, a) {}

//template <class ArgsTuple>
//double eval(const double &, const ArgsTuple &) const {
//    return get_value();
//}
//template <class ArgsTuple>
//double norm(const ArgsTuple &) const
//{
//    return std::abs(get_value());
//}
//int to_int() const {
//    if (!is_integer(get_value())) {
//        PIRANHA_THROW(value_error, "cannot convert double coefficient to integer");
//    }
//    return (int)get_value();
//}
//template <class ArgsTuple>
//std::complex<double_cf> ei(const ArgsTuple &) const;

// These tests are double_cf specific but are also using the correpsonding procedures and templates of NumericalContainer
BOOST_AUTO_TEST_CASE(construction_test)
{
    double_cf value01;
    BOOST_TEST(value01.get_value() == 0.0);

    //Argstuple is nowhere used why is it here??
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    double_cf value02(99.01, emptyTuple);
    BOOST_TEST(value02.get_value() == 99.01);

    mp_rational half("1/2");                             // I think only mp_rational and mp_integer are supported
    double_cf value03(half, emptyTuple);
    BOOST_TEST(value03.get_value() == 0.5);

    mp_integer one("1");                                   // I think only mp_rational and mp_integer are supported 
    double_cf value04(one, emptyTuple);
    BOOST_TEST(value04.get_value() == 1.0);

    // the value is 1, what is Psym and argTuple actually used for??. They are completely ignored
    Psym v1{ "v1" };
    double_cf value05(v1, 99, emptyTuple);
    BOOST_TEST(value05.get_value() == 1.0);

    // copy construction
    double_cf value06(value02);
    BOOST_TEST(value06 == value02);
    BOOST_TEST(value06.get_value() == 99.01);

    double_cf value07(string("-2.51"), emptyTuple); // does not work for char * !!, needs new definition
    BOOST_TEST(value07.get_value() == -2.51);
    // NUmericalContainer also has a constructor from string and argsTuple. Where is that used??? Not for double_cf.
}

BOOST_AUTO_TEST_CASE(eval_test)
{
    // identical to getValue(), just a wrapper but has parameters. NOne of the paramaeters is used!!

    //Argstuple is nowhere used why is it here??
    // what is the double in the parameter list good for???It is not used
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    BOOST_TEST_CONTEXT("return value construct from value double");  // has first parameter as template. How many to create??
    double_cf value02(99.01, emptyTuple);
    BOOST_TEST(value02.eval(100.0, emptyTuple) == 99.01);

    BOOST_TEST_CONTEXT("return construct from value rational"); // uses numerical container constructor selector
    mp_rational half("1/2");                             // I think only mp_rational and mp_integer are supported
    double_cf value03(half, emptyTuple);
    BOOST_TEST(value03.eval(100.0, emptyTuple) == 0.5);

    BOOST_TEST_CONTEXT("Construct from value mp integer"); // uses numerical container constructor selector
    mp_integer one("1");                                   // I think only mp_rational and mp_integer are supported 
    double_cf value04(one, emptyTuple);
    BOOST_TEST(value04.eval(100.0,emptyTuple) == 1.0);
}

BOOST_AUTO_TEST_CASE(norm_test)
{
       // argstuple is nowhere used
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    BOOST_TEST_CONTEXT("return value construct from value double");  // has first parameter as template. How many to create??
    double_cf value02(99.01, emptyTuple);
    BOOST_TEST(value02.norm(emptyTuple) == 99.01);

    double_cf value03(-99.01, emptyTuple);
    BOOST_TEST(value03.norm(emptyTuple) == 99.01);

    mp_rational half("1/2");                             // I think only mp_rational and mp_integer are supported
    double_cf value04(half, emptyTuple);
    BOOST_TEST(value04.norm(emptyTuple) == 0.5);

    mp_rational nhalf("-1/2");                             // I think only mp_rational and mp_integer are supported
    double_cf value05(half, emptyTuple);
    BOOST_TEST(value05.norm(emptyTuple) == 0.5);
}

BOOST_AUTO_TEST_CASE(to_integer_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value02(99.01, emptyTuple);
    BOOST_CHECK_THROW(value02.to_int(), value_error);

    double_cf value03(99.00, emptyTuple);
    BOOST_TEST(value03.to_int() == 99);
}

// test as an extension of std the double_cf complex values
//explicit complex() : ancestor() {}

//template <class T, class ArgsTuple>
//explicit complex(const T &x, const ArgsTuple &argsTuple) : ancestor(x, argsTuple) {}

//template <class ArgsTuple>
//explicit complex(const piranha::mp_rational &q, const ArgsTuple &argsTuple) : ancestor(q.to_double(), argsTuple) {}

//template <class ArgsTuple>
//explicit complex(const piranha::mp_integer &z, const ArgsTuple &argsTuple) : ancestor(z.to_double(), argsTuple) {}

//template <class ArgsTuple>
//explicit complex(const piranha::Psym &p, const int &n, const ArgsTuple &a) : ancestor(p, n, a) {}

//template <class ArgsTuple>
//complex<double> eval(const double &, const ArgsTuple &) const {
//    return get_value();
//}
//template <class ArgsTuple>
//double norm(const ArgsTuple &) const {
//    return abs(get_value());
//}

BOOST_AUTO_TEST_CASE(construction_test_complex)
{
    using complexCf = complex<double_cf>;
    BOOST_TEST_CONTEXT("Simple construction");
    complexCf value01;
    BOOST_TEST(value01.get_value().real() == 0.0);
    BOOST_TEST(value01.get_value().imag() == 0.0);

    //Argstuple is nowhere used why is it here??
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    BOOST_TEST_CONTEXT("Construct from value double");  // has first parameter as template. How many to create??
    complexCf value02(99.01, emptyTuple);
    BOOST_TEST(value02.get_value() == 99.01); // what does comparison precisely contain. comparison of complex with real!

    BOOST_TEST_CONTEXT("Construct from value rational"); // uses numerical container constructor selector
    mp_rational half("1/2");                             // I think only mp_rational and mp_integer are supported
    complexCf value03(half, emptyTuple);
    BOOST_TEST(value03.get_value() == 0.5);

    BOOST_TEST_CONTEXT("Construct from value mp integer"); // uses numerical container constructor selector
    mp_integer one("1");                                   // I think only mp_rational and mp_integer are supported 
    complexCf value04(one, emptyTuple);
    BOOST_TEST(value04.get_value() == 1.0);

    // the value is 1, what is Psym and argTuple actually used for??. They are completely ignored
    BOOST_TEST_CONTEXT("Construct from with Psymr");
    Psym v1{ "v1" };
    complexCf value05(v1, 99, emptyTuple);
    BOOST_TEST(value05.get_value() == 1.0);

    //double_cf value07(string("-2.51,3.01"), emptyTuple); // does not work for char * !!, needs new definition
    //BOOST_TEST(value07.get_value() == (-2.51,3.01));
    // NUmericalContainer also has a constructor from string and argsTuple. Where is that used??? Not for double_cf.

    // TODO: how to construct a complex from a complex literal value
}

