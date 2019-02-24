#define BOOST_TEST_MODULE Double coefficients
#include "boost/test/included/unit_test.hpp"
#include "boost/test/output_test_stream.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// contains also all test cases for the base class NumericalContainer
// Numerical Container is only used as base class for other classes and can not be instantaiated by itself.
//
// specific implementation in derived classes are tested in the corresponding tests.
// These tests for doubleCf serve as default tests for the NumericalContainer class
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include "piranha.h"

using namespace piranha;
using namespace std;
using boost::test_tools::output_test_stream;
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

    double_cf value07(string(" -2.51"), emptyTuple); // does not work for char * !!, needs new definition
    BOOST_TEST(value07.get_value() == -2.51);
}

BOOST_AUTO_TEST_CASE(eval_test)
{
    // identical to getValue(), just a wrapper but has parameters. None of the parameters is used!!

    //Argstuple is nowhere used why is it here??
    // what is the double in the parameter list good for???It is not used
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    double_cf value02(99.01, emptyTuple);
    BOOST_TEST(value02.eval(100.0, emptyTuple) == 99.01);

    mp_rational half("1/2");                             // I think only mp_rational and mp_integer are supported
    double_cf value03(half, emptyTuple);
    BOOST_TEST(value03.eval(100.0, emptyTuple) == 0.5);

    mp_integer one("1");                                   // I think only mp_rational and mp_integer are supported 
    double_cf value04(one, emptyTuple);
    BOOST_TEST(value04.eval(100.0,emptyTuple) == 1.0);
}

BOOST_AUTO_TEST_CASE(norm_test)
{
       // argstuple is nowhere used
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
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
    complexCf value01;
    BOOST_TEST(value01.get_value().real() == 0.0);
    BOOST_TEST(value01.get_value().imag() == 0.0);

    //Argstuple is nowhere used why is it here??
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    complexCf value02(99.01, emptyTuple);
    BOOST_TEST(value02.get_value() == 99.01); // what does comparison precisely contain. comparison of complex with real!

    mp_rational half("1/2");                             // I think only mp_rational and mp_integer are supported
    complexCf value03(half, emptyTuple);
    BOOST_TEST(value03.get_value() == 0.5);

    mp_integer one("1");                                   // I think only mp_rational and mp_integer are supported 
    complexCf value04(one, emptyTuple);
    BOOST_TEST(value04.get_value() == 1.0);

    // the value is 1, what is Psym and argTuple actually used for??. They are completely ignored
    Psym v1{ "v1" };
    complexCf value05(v1, 99, emptyTuple);
    BOOST_TEST(value05.get_value() == 1.0);

    //double_cf value07(string("-2.51,3.01"), emptyTuple); // does not work for char * !!, needs new definition
    //BOOST_TEST(value07.get_value() == (-2.51,3.01));
    complex<double> com(1.0, -1.0);
    //complexCf value06(com);   // TODO: can't create complexCf from complex<double>???

    // TODO: how to construct a complex from a complex literal value? not implemented
}




///// Print in plain mode.
//template <class ArgsTuple>
//void printPlain(std::ostream &outStream, const ArgsTuple &) const {
//    outStream << boost::lexical_cast<std::string>(m_value);
//}
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// most of the subsequent tests are actually tests for NumaericalContainer. Some have specific implementations
// in derived classes. If so we have tests for the derived class
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(printPlain_test)
{
    // has specific implementation in cf_series
    
    //Argstuple is nowhere used why is it here??
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    output_test_stream output;
    double_cf value02(-99.0, emptyTuple);
    value02.printPlain(output, emptyTuple); // argsTuple nor used?
    BOOST_TEST(!output.is_empty(false));
    BOOST_TEST(output.is_equal("-99")); // not really a good test. How can we control better the output format ????
}

///// Print in pretty mode. Equivalent to print_plain.
//template <class ArgsTuple>
//void printPretty(std::ostream &outStream, const ArgsTuple &) const {
//    numerical_container_print_pretty_selector<T>::run(m_value, outStream);
//}
BOOST_AUTO_TEST_CASE(printPretty_test)
{
    // has specific implementation in cf_series. How about complex and e.g. mpq_cf as paramaeter?

//Argstuple is nowhere used why is it here??
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    output_test_stream output;
    double_cf value02(-99.0, emptyTuple); // there is not enough control over the output format which makes double tests fail
    value02.printPretty(output, emptyTuple); // argsTuple nor used?
    BOOST_TEST(!output.is_empty(false));
    BOOST_TEST(output.is_equal("-99")); // not a good test, but no real control over output format.
}


//
///// Print in Tex mode.
//template <class ArgsTuple>
//void printTex(std::ostream &outStream, const ArgsTuple &) const
//{
//    numerical_container_print_tex_selector<T>::run(m_value, outStream);
//}
BOOST_AUTO_TEST_CASE(printTex_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    output_test_stream output;
    double_cf value02(-99.0, emptyTuple); // there is not enough control over the output format which makes double tests fail
    value02.printTex(output, emptyTuple); // argsTuple nor used?
    BOOST_TEST(!output.is_empty(false));
    BOOST_TEST(output.is_equal(" -99 ")); // not a good test, but no real control over output format. Observe the blanks!
}
///// Swap content using std::swap.
//Derived &swap(Derived &dc)
//{
//    std::swap(m_value, dc.m_value);
//    return *derived_cast;
//}
BOOST_AUTO_TEST_CASE(swap_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    double_cf value00(1.5, emptyTuple); // why do we not have a constructor for just double value??
    double_cf value01(-9.1, emptyTuple);
    value00.swap(value01);
    BOOST_TEST(value00.get_value() == -9.1);
    BOOST_TEST(value01.get_value() == 1.5);
}
///// Pad right.
//template <class ArgsTuple>
//void padRight(const ArgsTuple &) {}
BOOST_AUTO_TEST_CASE(padRight_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    // does nothing, has specific implementation for series_cf
    double_cf value;
    value.padRight(emptyTuple); 
}

///// Apply layout.
//template <class Layout, class ArgsTuple>
//void applyLayout(const Layout &, const ArgsTuple &) {}
//
BOOST_AUTO_TEST_CASE(applyLayout_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value;

    LayoutTuple<boost::tuple<VectorPsym>>::Type layoutTuple;
    // no element at all
    LayoutElement nullPair = std::make_pair(false, 3); // LayoutElement is defined in named series!!!!
    Layout layout;
    layoutTuple.get<0>() = Layout();


    value.applyLayout(layout, emptyTuple); // does nothing what is it for? Do we need layout on this? series_cf maybe?
}

///// Test if trimming is possible.
//template <class TrimFlags>
//void trimTest(TrimFlags &) const {}

BOOST_AUTO_TEST_CASE(trimTest_test)
{
    using TrimFlagsType = NTuple<std::vector<bool>, 1>::Type;
    std::vector<bool> flags(1, true);
    TrimFlagsType trimFlags;
    trimFlags.get<0>() = flags;

    double_cf value;
    value.trimTest(trimFlags); // does nothing i.e. should return the flags as they are. 
                               // shouldn't it throw (attempt to trim but nothing to) or declare everything trimmable i.e. return false 
                                // does series_cf use it???
    BOOST_TEST(trimFlags.get<0>()[0] == true);
}

///// Trim.
//template <class TrimFlags, class ArgsTuple>
//Derived trim(const TrimFlags &, const ArgsTuple &) const
//{
//    return *derived_const_cast;
//}
BOOST_AUTO_TEST_CASE(trim_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    using TrimFlagsType = NTuple<std::vector<bool>, 1>::Type;
    std::vector<bool> flags(1, true);
    TrimFlagsType trimFlags;
    trimFlags.get<0>() = flags;

    double_cf value(-1.1, emptyTuple);
    BOOST_TEST(value.trim(trimFlags, emptyTuple) == value); // trim does nothing. Should it do something? when is it actually really usefull. Only if _cf
                                                            // has key that can be trimmend. series_cf.
}


///// Split.
//template <class Series, class ArgsTuple>
//void split(std::vector<std::vector<Series> > &, const int &, const ArgsTuple &) const
//{
//    PIRANHA_ASSERT(false);
//}
BOOST_AUTO_TEST_CASE(A_split_test)
{
    // needs serie. Why is defined on this level???    
    double_cf value00;
    //BOOST_CHECK_THROW(value00.split(), assertion_error); // not usefull at all. only for series_cf. Why is it on this level????
    BOOST_TEST_FAIL("Not implemented. Unclear if test belongs here");
}

///// Number of atoms. Returns 1.
//std::size_t atoms() const {
//    return 1;
//}

BOOST_AUTO_TEST_CASE(atoms_test)
{
    double_cf value00;
    BOOST_TEST(value00.atoms() == 1); // where is that actually used? Here it is always 1. Why? Shouldn't it implemented by the actual implementation??
}
///// Test for ignorability.
///**
// * Returns true if norm() is less than settings::get_numerical_zero().
// */
//template <class ArgsTuple>
//bool isIgnorable(const ArgsTuple &a) const
//{
//    return (derived_const_cast->norm(a) < settings::get_numerical_zero());
//}
BOOST_AUTO_TEST_CASE(isIgnorable_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    
    double_cf value00(settings::get_numerical_zero() * 2, emptyTuple); // this depends on settings!!!!!
    BOOST_TEST(value00.isIgnorable(emptyTuple) == false);

    double_cf value01(settings::get_numerical_zero() / 2, emptyTuple);// what should this be. This depends heavily on the chosen limits.
                                                                      // DOes it make sense when we can deal mpr high precission numerics????
    BOOST_TEST(value01.isIgnorable(emptyTuple) == true);

    double_cf value02(0, emptyTuple);
    BOOST_TEST(value02.isIgnorable(emptyTuple) == true);
}

///// Insertability test. Returns true.
//template <class ArgsTuple>
//bool isInsertable(const ArgsTuple &) const
//{
//    return true;
//}
BOOST_AUTO_TEST_CASE(isInsertable_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    // it is always insertable except for ??? series_cf??
    double_cf value00;
    BOOST_TEST(value00.isInsertable(emptyTuple) == true);
}

///// Padding test. Returns false.
//template <class ArgsTuple>
//bool needsPadding(const ArgsTuple &) const
//{
//    return false;
//}
BOOST_AUTO_TEST_CASE(needsPadding_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    double_cf value00;
    BOOST_TEST(value00.needsPadding(emptyTuple) == false); // always false. When is it true???
}

//template <class U>
//bool operator==(const U &x) const
//{
//    return numerical_container_equality_selector<Derived, U>::run(*derived_const_cast, x);
//}
BOOST_AUTO_TEST_CASE(equality_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    
    double_cf value00(9.99, emptyTuple);
    double_cf value01(1.99, emptyTuple);
    double_cf value02(9.99, emptyTuple);
    BOOST_TEST(value00 == 9.99);   // where and how are these comparisons defined??
    BOOST_TEST(!(value00 == 0.0)); 
    BOOST_TEST(value00 == value02);
    BOOST_TEST(!(value00 == value01));
}
//template <class U>
//bool operator!=(const U &x) const
//{
//    return !operator==(x);
//}
BOOST_AUTO_TEST_CASE(nonequality_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    double_cf value00(9.99, emptyTuple);
    double_cf value02(-9.99, emptyTuple);
    BOOST_TEST(value00 != -99.0);
    BOOST_TEST(value00 != value02);
}
//// Math.
//template <class ArgsTuple>
//void invertSign(const ArgsTuple &)
//{
//    m_value *= -1;
//}

BOOST_AUTO_TEST_CASE(invertSign_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(1.0, emptyTuple);
    value00.invertSign(emptyTuple);
    BOOST_TEST(value00.get_value() == -1.0);
}
//template <class U, class ArgsTuple>
//Derived &add(const U &x, const ArgsTuple &)
//{
//    return numerical_container_add_selector<U>::run(*derived_cast, x);
//}
BOOST_AUTO_TEST_CASE(add_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(1.0, emptyTuple);
    
    double_cf value01 = value00.add(2.0, emptyTuple); // uses copy constructor
    BOOST_TEST(value01.get_value() == 3.0);
    double_cf value02;
    value02 = value00.add(3.0, emptyTuple); // uses assignment operator. Where is that defined????
    BOOST_TEST(value02.get_value() == 6.0); // add is accumulative on value00
}
//template <class U, class ArgsTuple>
//Derived &subtract(const U &x, const ArgsTuple &)
//{
//    return numerical_container_subtract_selector<U>::run(*derived_cast, x);
//}
BOOST_AUTO_TEST_CASE(subtract_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(1.0, emptyTuple);

    double_cf value01 = value00.subtract(2.0, emptyTuple); // uses copy constructor
    BOOST_TEST(value01.get_value() == -1.0);
    double_cf value02;
    value02 = value00.subtract(3.0, emptyTuple); // uses assignment operator. Where is that defined????
    BOOST_TEST(value02.get_value() == -4.0);  // sub is accumulative on value00

    // using "-" directly is not defined. What is the argutple good for???
    //double_cf value03;
    //value03 = value00 - value02;
}
//template <class U, class ArgsTuple>
//Derived &multBy(const U &x, const ArgsTuple &)
//{
//    return numerical_container_multiply_selector<U>::run(*derived_cast, x);
//}
BOOST_AUTO_TEST_CASE(multBy_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(2.0, emptyTuple);

    value00.multBy(-3.0, emptyTuple);

    BOOST_TEST(value00.get_value() == -6.0);

    double_cf value01(2.0, emptyTuple);
    double_cf value02(-2.0, emptyTuple);

    value00 = value01.multBy(value02, emptyTuple);
    BOOST_TEST(value01.get_value() == -4.0);
    BOOST_TEST(value00.get_value() == -4.0);
}

//template <class U, class ArgsTuple>
//Derived &divideBy(const U &x, const ArgsTuple &)
//{
//    return numerical_container_divide_selector<U>::run(*derived_cast, x);
//}
BOOST_AUTO_TEST_CASE(divideBy_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(2.0, emptyTuple);

    value00.divideBy(-3.0, emptyTuple);

    BOOST_TEST(value00.get_value() == -2.0/3.0);

    double_cf value01(2.0, emptyTuple);
    double_cf value02(-2.0, emptyTuple);

    value00 = value01.divideBy(value02, emptyTuple);
    BOOST_TEST(value01.get_value() == -1.0);
    BOOST_TEST(value00.get_value() == -1.0);

}
//// Multiply and add.
//template <class Derived2, class ArgsTuple>
//void addmul(const Derived &x1, const Derived2 &x2, const ArgsTuple &)
//{
//    multiply_accumulate(m_value, x1.m_value, x2.get_value());
//}
// calculates a*b + c (speed)
BOOST_AUTO_TEST_CASE(addmul_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(2.0, emptyTuple); //c
    double_cf value01(2.0, emptyTuple); //a
    double_cf value02(-3.0, emptyTuple);//b
    value00.addmul(value01, value02, emptyTuple);
    BOOST_TEST(value00.get_value() == -4.0);
}

//template <class Series, class PosTuple, class ArgsTuple>
//Series partial(const PosTuple &, const ArgsTuple &) const
//{
//    return Series();
//}
// why is the definition on this level
BOOST_AUTO_TEST_CASE(A_partial_test)
{
    // requires a series but here anything can serve. There is no test on it
    //VectorPsym noSymbol;
    //boost::tuple<VectorPsym> emptyTuple(noSymbol);
    //double_cf value00(2.0, emptyTuple); //c

    //double_cf value01 = value01.partial<double_cf,int,int>(2, 2, 2);
    BOOST_FAIL("Not implemented. Unclear if test belongs here");
}

//template <class ArgsTuple>
//Derived pow(const double &x, const ArgsTuple &argsTuple) const
//{
//    return generic_pow(x, argsTuple);
//}
BOOST_AUTO_TEST_CASE(powdouble_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(2.0, emptyTuple);
    double_cf value01;
    value01 = value00.pow(2.0, emptyTuple); // other than add/mul/sub power only returns the operated value
    BOOST_TEST(value01.get_value() == 4.0);
}

//template <class ArgsTuple>
//Derived pow(const mp_rational &q, const ArgsTuple &argsTuple) const
//{
//    return generic_pow(q, argsTuple);
//}
BOOST_AUTO_TEST_CASE(powrational_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(2.0, emptyTuple);
    double_cf value01;
    value01 = value00.pow(mp_rational("2"), emptyTuple);
    BOOST_TEST(value01.get_value() == 4.0);
}
//template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
//RetSeries sub(const PosTuple &, SubCaches &, const ArgsTuple &argsTuple) const
//{
//    return RetSeries::baseSeriesFromCf(*derived_const_cast, argsTuple);
//}
BOOST_AUTO_TEST_CASE(A_substitute_test)
{
    // requires series. Why is on this level
    BOOST_FAIL("Not implemented. Unclear if test belongs here");
}
//template <class RetSeries, class PosTuple, class SubCaches, class ArgsTuple>
//RetSeries eiSubstitute(const PosTuple &p, SubCaches &s, const ArgsTuple &a) const
//{
//    return sub<RetSeries>(p, s, a);
//}
BOOST_AUTO_TEST_CASE(A_eisubstitute_test)
{
    // requires series. Why is on this level
    BOOST_FAIL("Not implemented. Unclear if test belongs here");
}
///// Bessel function of the first kind.
///**
// * Will use piranha::besselJ internally.
// */
//template <class ArgsTuple>
//Derived besselJ(const int &n, const ArgsTuple &argsTuple) const
//{
//    return Derived(piranha::besselJ(n, m_value), argsTuple);
//}
BOOST_AUTO_TEST_CASE(besselJ_test)
{
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(3.0, emptyTuple);

    double_cf value01;
    value01 = value00.besselJ(1, value00);
    BOOST_TEST(value01.get_value() == 0.33905895852593643);
    value01 = value00.besselJ(2, value00);
    BOOST_TEST(value01.get_value() == 0.48609126058589103);
}
///// Get value.
//const T &get_value() const
//{
//    return m_value;
//}
BOOST_AUTO_TEST_CASE(getValue_test)
{   
    // see below. we don't use set_value!
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);

    double_cf value00(-9.99, emptyTuple);
    BOOST_TEST(value00.get_value() == -9.99);
}

///// Set value.
//template <class U>
//void set_value(const U &value)
//{
//    m_value = value;
//}
BOOST_AUTO_TEST_CASE(setValue_test) // why do we have this but not an assignment operator????
{    
    double_cf value00;
    value00.set_value(-1);
    BOOST_TEST(value00.get_value() == -1);
}

///// get as string
//
//inline std::string toString() const
//{
//    return boost::lexical_cast<std::string>(m_value);
//}
BOOST_AUTO_TEST_CASE(toString_test)
{
    // not properly testable. depends on the output string formatting which is not properly defined!!
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(2.0, emptyTuple);

    BOOST_TEST(value00.toString() == "2");
}

//inline std::size_t printLength() const
//{
//    return toString().length();
//}
BOOST_AUTO_TEST_CASE(printlength_test)
{
    // very fuzzy. It all depends on the formatting. there is no real defined length. depends
    // on the output formatting which is not properly defined. This is an artifically adjusetd test
    VectorPsym noSymbol;
    boost::tuple<VectorPsym> emptyTuple(noSymbol);
    double_cf value00(2.0, emptyTuple);
    BOOST_TEST(value00.toString().length() == value00.printLength()); // this should always be true
}


