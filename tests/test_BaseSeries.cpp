#define BOOST_TEST_MODULE BaseSeries Test
#include "boost/test/included/unit_test.hpp"
#include "boost/test/tools/output_test_stream.hpp"


#include "piranha.h"

using namespace std;
using namespace piranha;
using boost::test_tools::output_test_stream;

//
// These are tests for BaseSeries. The class BaseSeries is not intended to be used on it's own but we can test
// it independently. It is supposed to be inhereted e.g. by named series etc to provide the lowest level of series/polynomial
// manipulation, but we can test this functionality without all the additional clases to provide a full "blown" e.g. polynomial.
// For this purpose we use a simple version of it to minimize the effort
//
// For this purpose we create an intermediate class, which is practically a "polynomial" without having any real polynomial features
// but to be able to test the primitive functions without all the overhead comming from a "real" polynomial.
// This is similar to what we do for VectorKey tests.


// argTuples are introduced using Psym but their connection to names is not actually being used! That is bussines of named series
// Does that smell like bad design??

namespace {

    template< class Term>
    class TSeries : public BaseSeries < Term, '!', std::allocator<char>, TSeries<Term> > 
    {
    public: 
        TSeries() = default;
        TSeries(TSeries const&) = default;
        TSeries(TSeries&&) = default;

        TSeries& operator=(TSeries const&) = default;
        TSeries& operator=(TSeries&&) = default;

        virtual ~TSeries() = default;
    };

    using STerm = Monomial<double_cf, ExpoVector<int, 0>, '|', std::allocator<char>>; // a ploynomial like term with double coefficients
    using SSeries = TSeries<STerm>;
}


BOOST_AUTO_TEST_CASE(defaultconstruction_test)
{
    SSeries series;
    BOOST_TEST(series.length() == 0);
    BOOST_TEST(series.echelonLevel == 0);
    BOOST_TEST(series.separator == '!');
    BOOST_TEST(series.empty() == true);
    BOOST_TEST(series.isSingleCf() == false);
    BOOST_TEST(series.atoms() == 0);
}

//template <class Key, class ArgsTuple>
//static Derived baseSeriesFromKey(const Key&, const ArgsTuple&);
BOOST_AUTO_TEST_CASE(BaseSeriesFromKeyTest)
{
    Psym p1("p1");
    Psym p2("p2");
    Psym p3("p3");
    VectorPsym polySym = { p1, p2, p3 };
    boost::tuple<VectorPsym> polyOnlyArgs(polySym);

    SSeries::TermType::KeyType key;
    key.resize(3);
    key[0] = 1;
    key[1] = 2;
    key[2] = 3;
    SSeries bseries01 = SSeries::baseSeriesFromKey(key, polyOnlyArgs);
    BOOST_TEST(bseries01.length() == 1);
    BOOST_TEST(bseries01.echelonLevel == 0);
    BOOST_TEST(bseries01.empty() == false);
    BOOST_TEST(bseries01.isSingleCf() == false);
    BOOST_TEST(bseries01.atoms() == 2);

    SSeries::TermType term = *(bseries01.begin());
    BOOST_TEST(term.key[0] == 1);
    BOOST_TEST(term.key[1] == 2);
    BOOST_TEST(term.key[2] == 3);
    BOOST_TEST(term.cf == 1);

	// no key coeficients which means it is only a coefficient
    key.resize(3);
    key[0] = 0;
    key[1] = 0;
    key[2] = 0;
    SSeries bseries02 = SSeries::baseSeriesFromKey(key, polyOnlyArgs);
    BOOST_TEST(bseries02.length() == 1);
    BOOST_TEST(bseries02.echelonLevel == 0);
    BOOST_TEST(bseries02.empty() == false);
    BOOST_TEST(bseries02.isSingleCf() == true);
    BOOST_TEST(bseries02.atoms() == 2);

    term = *(bseries02.begin());
    BOOST_TEST(term.key[0] == 0);
    BOOST_TEST(term.key[1] == 0);
    BOOST_TEST(term.key[2] == 0);
    BOOST_TEST(term.cf == 1);


    // how to test the element in the series???
}

//
//template <class Cf, class ArgsTuple>
//static Derived baseSeriesFromCf(const Cf&, const ArgsTuple&);

 BOOST_AUTO_TEST_CASE(BaseSeriesFromCfTest)
 {
     Psym p1("p1");
     Psym p2("p2");
     Psym p3("p3");
     VectorPsym polySym = { p1, p2, p3 };
     boost::tuple<VectorPsym> polyOnlyArgs(polySym);

     // *cf can be all kinds of types, SHould be implemented
     SSeries bseries03 = SSeries::baseSeriesFromCf(double_cf(2.0, polyOnlyArgs), polyOnlyArgs);

     BOOST_TEST(bseries03.length() == 1);
     BOOST_TEST(bseries03.echelonLevel == 0);
     BOOST_TEST(bseries03.empty() == false);
     BOOST_TEST(bseries03.isSingleCf() == true);
     BOOST_TEST(bseries03.atoms() == 2);

     SSeries::TermType term = *(bseries03.begin());
     BOOST_TEST(term.key[0] == 0);
     BOOST_TEST(term.key[1] == 0);
     BOOST_TEST(term.key[2] == 0);
     BOOST_TEST(term.cf == 2.0);

     // how to test the element in the series???
 }
//template <class Number, class ArgsTuple>
//static Derived baseSeriesFromNumber(const Number&, const ArgsTuple&);
 BOOST_AUTO_TEST_CASE(BaseSeriesFromNumber)
 {
     Psym p1("p1");
     Psym p2("p2");
     Psym p3("p3");
     VectorPsym polySym = { p1, p2, p3 };
     boost::tuple<VectorPsym> polyOnlyArgs(polySym);

     SSeries bseries04 = SSeries::baseSeriesFromNumber(2.0, polyOnlyArgs);

     BOOST_TEST(bseries04.length() == 1);
     BOOST_TEST(bseries04.echelonLevel == 0);
     BOOST_TEST(bseries04.empty() == false);
     BOOST_TEST(bseries04.isSingleCf() == true);
     BOOST_TEST(bseries04.atoms() == 2);

     // how to test the element in the series???
     SSeries::TermType term = *(bseries04.begin());
     BOOST_TEST(term.key[0] == 0);
     BOOST_TEST(term.key[1] == 0);
     BOOST_TEST(term.key[2] == 0);
     BOOST_TEST(term.cf == 2.0);
 }

//template <class ArgsTuple>
//void baseConstructFromPsym(const Psym&, const int, const ArgsTuple&);
 // Not sure that this is used anywhere explicitly
 // there is reference to it in named series via "constructFromPsym" which doesn't seem to be used
 // anywhere as well as in constructor polynomial_cf which also doesn;t seem to be in use
 // i.e. purpose ?? and even is it correct (which requires a purpose)
 BOOST_AUTO_TEST_CASE(baseConstructFromPsym)
 {
     Psym p1("p1");
     Psym p2("p2");
     Psym p3("p3");
     VectorPsym polySym = { p1, p2, p3 };
     boost::tuple<VectorPsym> polyOnlyArgs(polySym);

     SSeries bseries05;
     // has to be the echelon level of the key (=0) otherwise nothing is done (empty series)
     // the args tuple has to have only one element adn this element has to agree with the psym
     // STill: what is that good for?
     // this creates a series with 1 term of value 1
     bseries05.baseConstructFromPsym(p1, 1, polyOnlyArgs); 

     BOOST_TEST(bseries05.length() == 1);
     BOOST_TEST(bseries05.echelonLevel == 0);
     BOOST_TEST(bseries05.empty() == false);
     BOOST_TEST(bseries05.isSingleCf() == true); // key is unity

     // how to test the element in the series???
     SSeries::TermType term = *(bseries05.begin());
     BOOST_TEST(term.key[0] == 0);
     BOOST_TEST(term.key[1] == 0);
     BOOST_TEST(term.key[2] == 0);
     BOOST_TEST(term.cf == 1.0); // ???? why doe this not fail

     SSeries bseries06;

     VectorPsym polySym1 = { p1 };
     boost::tuple<VectorPsym> polyOnlyArgs1(polySym1);

     // should throw assertion because p2 is not in polyOnlyArgs1
     BOOST_CHECK_THROW(bseries06.baseConstructFromPsym(p2, 0, polyOnlyArgs1), assertion_error); // echelon level is 0

     SSeries bseries07;
     // should throw assertion because is not in polyOnlyArgs  has more than 1 element
     BOOST_CHECK_THROW(bseries07.baseConstructFromPsym(p1, 0, polyOnlyArgs), assertion_error); // echelon level is 0

     SSeries bseries08;
     // should go through because p1 is the only element in polyOnlyArgs1 
     BOOST_CHECK_NO_THROW(bseries08.baseConstructFromPsym(p1, 0, polyOnlyArgs1)); // echelon level is 0

     BOOST_TEST(bseries08.length() == 1);
     BOOST_TEST(bseries08.echelonLevel == 0);
     BOOST_TEST(bseries08.empty() == false);
     BOOST_TEST(bseries08.isSingleCf() == false); // key is not unity

     term = *(bseries08.begin());
     BOOST_TEST(term.key.size() == 1);
     BOOST_TEST(term.key[0] == 1);
     BOOST_TEST(term.cf == 1.0);
 }


 //void eraseTerm(const const_iterator&);
 // only by iterator. By term would be helpful too.
 BOOST_AUTO_TEST_CASE(eraseTerm)
 {
	 Psym p1("p1");
	 Psym p2("p2");
	 Psym p3("p3");
	 VectorPsym polySym = { p1, p2, p3 };
	 boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	 SSeries::TermType::KeyType key01;
	 key01.resize(3);

	 // create a vector of 5 terms
	 std::vector<SSeries::TermType> sourceVector;
	 for (int i = 0; i < 5; ++i)
	 {
		 // we reuse the key
		 key01[0] = 1 * i;
		 key01[1] = 2 * i;
		 key01[2] = -2 * i;

		 SSeries::TermType::CfType value01(1.0 * i, polyOnlyArgs);
		 sourceVector.push_back(SSeries::TermType(value01, key01));
	 }

	 BOOST_TEST(sourceVector.size() == 5);  // just to verify we actually have 5 terms

	 SSeries bseries01; // the baseseries we populate
	 for (int i = 0; i < 5; ++i)
	 {
		 BOOST_CHECK_NO_THROW(bseries01.insert(sourceVector[i], polyOnlyArgs));
	 }

	 BOOST_TEST(bseries01.length() == 4); // all elements are different should have length of vector-1.
										  // index 0 constructs a null term(!) and is discarded	

	 // index 0 term shall not be in the series (a null term)
	 BOOST_TEST((bseries01.findTerm(sourceVector[0]) == bseries01.end()));

	 for (int i = 1; i < 5; ++i)
	 {
		 bseries01.eraseTerm(bseries01.findTerm(sourceVector[i]));
		 BOOST_TEST((bseries01.findTerm(sourceVector[i]) == bseries01.end()));
		 BOOST_TEST((bseries01.length() == 4 - i));
	 }
 }

 //void clearTerms();
 BOOST_AUTO_TEST_CASE(clearTerms)
 {
	 Psym p1("p1");
	 Psym p2("p2");
	 Psym p3("p3");
	 VectorPsym polySym = { p1, p2, p3 };
	 boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	 SSeries::TermType::KeyType key01;
	 key01.resize(3);

	 // create a vector of 5 terms
	 SSeries bseries01; // the baseseries we populate
	 for (int i = 0; i < 5; ++i)
	 {
		 // we reuse the key
		 key01[0] = 1 * i;
		 key01[1] = 2 * i;
		 key01[2] = -2 * i;

		 SSeries::TermType::CfType value01(1.0 * i, polyOnlyArgs);
		 bseries01.insert(SSeries::TermType(value01, key01), polyOnlyArgs);
	 }

	 BOOST_TEST(bseries01.length() == 4);

	 bseries01.clearTerms();
	 BOOST_TEST(bseries01.length() == 0);
 }

 //template <bool, bool, class Term2, class ArgsTuple>
 //void insert(const Term2&, const ArgsTuple&);
 // insert with canonical check(template paramter 1)  and/or sign (template parameter 2)
 // dependent on parameters differetn operations are executed
 BOOST_AUTO_TEST_CASE(insertTest)
 {
     Psym p1("p1");
     Psym p2("p2");
     Psym p3("p3");
     VectorPsym polySym = { p1, p2, p3 };
     boost::tuple<VectorPsym> polyOnlyArgs(polySym);

     SSeries bseries01;

     SSeries::TermType::KeyType key01("1;2", polyOnlyArgs);
     SSeries::TermType::CfType value01(-1.0, polyOnlyArgs);

     SSeries::TermType term01(value01, key01); // how about something like term01(1.0,"1,2,3") etc. ???

     // bseries01 gets extended to number of key elements as they are described by polyOnlyArgs
     BOOST_CHECK_NO_THROW((bseries01.insert<false, false>(term01, polyOnlyArgs))); // no canonical check, no sign which means insertion of neagted term!

     BOOST_TEST(bseries01.length() = 1);
     BOOST_TEST(bseries01.isSingleCf() == false);

     SSeries::TermType term11 = *(bseries01.begin());
     BOOST_TEST(!(term11 == term01)); // different number of arguments
     // but
     BOOST_TEST(term11.key.size() == 3);
     BOOST_TEST(term11.key[0] == 1);
     BOOST_TEST(term11.key[1] == 2);
     BOOST_TEST(term11.key[2] == 0);
     BOOST_TEST(term11.cf == 1.0);


     SSeries bseries02;

     SSeries::TermType::KeyType key2("1;2", polyOnlyArgs);
     SSeries::TermType::CfType value02(-1.0, polyOnlyArgs);

     SSeries::TermType term02(value02, key2); // how about something like term01(1.0,"1,2,3") etc. ???

     // bseries01 gets extended to number of key elements as they are described by polyOnlyArgs
     BOOST_CHECK_NO_THROW((bseries02.insert<true, false>(term02, polyOnlyArgs))); // canonical check, no sign (i.e. insert inverted term) 
     // canonical check makes no difference we are dealing here with expoVector as key

     BOOST_TEST(bseries02.length() = 1);
     BOOST_TEST(bseries02.isSingleCf() == false);

     SSeries::TermType term12 = *(bseries02.begin());
     BOOST_TEST(term12.key.size() == 3);
     BOOST_TEST(term12.key[0] == 1);
     BOOST_TEST(term12.key[1] == 2);
     BOOST_TEST(term12.key[2] == 0); // this element did not exists originally. 
     BOOST_TEST(term12.cf == 1.0);


     SSeries bseries03;

     SSeries::TermType::KeyType key3("1;2", polyOnlyArgs);
     SSeries::TermType::CfType value03(-1.0, polyOnlyArgs);

     SSeries::TermType term03(value03, key3); // how about something like term01(1.0,"1,2,3") etc. ???

     // bseries01 gets extended to number of key elements as they are described by polyOnlyArgs
     BOOST_CHECK_NO_THROW((bseries03.insert<true, true>(term03, polyOnlyArgs))); // canonical check, with sign i.e. sddition!
     // canonical check makes no difference we are dealing here with expoVector as key

     BOOST_TEST(bseries03.length() = 1);
     BOOST_TEST(bseries03.isSingleCf() == false);

     SSeries::TermType term13 = *(bseries03.begin());
     BOOST_TEST(term13.key.size() == 3);
     BOOST_TEST(term13.key[0] == 1);
     BOOST_TEST(term13.key[1] == 2);
     BOOST_TEST(term13.key[2] == 0);
     BOOST_TEST(term13.cf == -1.0);


     // insert the same key twice
     // no sign (= false) i.e. subtract
     SSeries bseries04;

     SSeries::TermType::KeyType key4("1;2;3", polyOnlyArgs);
     SSeries::TermType::CfType value04(-1.0, polyOnlyArgs);

     SSeries::TermType term04(value04, key4); // how about something like term01(1.0,"1,2,3") etc. ???

     SSeries::TermType::CfType value24(-5.0, polyOnlyArgs);
     SSeries::TermType term24(value24, key4);
     // bseries01 gets extended to number of key elements as they are described by polyOnlyArgs
     BOOST_CHECK_NO_THROW((bseries04.insert<true, false>(term04, polyOnlyArgs))); // canonical check, without sign i.e. subtraction
     BOOST_CHECK_NO_THROW((bseries04.insert<true, false>(term24, polyOnlyArgs))); // canonical check, without sign i.e. subtraction

     // canonical check makes no difference we are dealing here with expoVector as key

     BOOST_TEST(bseries04.length() = 1);
     BOOST_TEST(bseries04.isSingleCf() == false);

     SSeries::TermType term14 = *(bseries04.begin());
     BOOST_TEST(term14.key.size() == 3);
     BOOST_TEST(term14.key[0] == 1);
     BOOST_TEST(term14.key[1] == 2);
     BOOST_TEST(term14.key[2] == 3);
     BOOST_TEST(term14.cf == 6.0);

 }

 //template <class Term2, class ArgsTuple>
 //void insert(const Term2&, const ArgsTuple&);
 // the standard insert function.
 // insert a term with canonicity check and sign (corresponds to just add it)
 BOOST_AUTO_TEST_CASE(insertSTdTest)
 {

	 // insert a single new term
	 Psym p1("p1");
	 Psym p2("p2");
	 Psym p3("p3");
	 VectorPsym polySym = { p1, p2, p3 };
	 boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	 SSeries bseries01;

	 SSeries::TermType::KeyType key01("1;2", polyOnlyArgs);
	 SSeries::TermType::CfType value01(-1.0, polyOnlyArgs);


	 SSeries::TermType term01(value01, key01); // 

	 // bseries01 gets extended to number of key elements as they are described by polyOnlyArgs
	 BOOST_CHECK_NO_THROW((bseries01.insert(term01, polyOnlyArgs))); // no canonical check, no sign which means insertion of neagted term!

	 BOOST_TEST(bseries01.length() = 1);
	 BOOST_TEST(bseries01.isSingleCf() == false);

	 SSeries::TermType term11 = *(bseries01.begin());
	 BOOST_TEST(!(term11 == term01)); // different number of arguments
	 // but
	 BOOST_TEST(term11.key.size() == 3);
	 BOOST_TEST(term11.key[0] == 1);
	 BOOST_TEST(term11.key[1] == 2);
	 BOOST_TEST(term11.key[2] == 0);
	 BOOST_TEST(term11.cf == -1.0);


	 //insert two terms with the same key
	 SSeries bseries04;

	 SSeries::TermType::KeyType key4("1;2;3", polyOnlyArgs);
	 SSeries::TermType::CfType value04(-1.0, polyOnlyArgs);

	 SSeries::TermType term04(value04, key4); // how about something like term01(1.0,"1,2,3") etc. ???

	 SSeries::TermType::CfType value24(-5.0, polyOnlyArgs);
	 SSeries::TermType term24(value24, key4);
	 // bseries01 gets extended to number of key elements as they are described by polyOnlyArgs
	 BOOST_CHECK_NO_THROW((bseries04.insert(term04, polyOnlyArgs))); // canonical check, and sign i.e. addition
	 BOOST_CHECK_NO_THROW((bseries04.insert(term24, polyOnlyArgs))); // canonical check, and sign i.e. addition

	 // canonical check makes no difference we are dealing here with expoVector as key
	 // should create a test case for TrigKey

	 BOOST_TEST(bseries04.length() = 1);
	 BOOST_TEST(bseries04.isSingleCf() == false);

	 SSeries::TermType term14 = *(bseries04.begin());
	 BOOST_TEST(term14.key.size() == 3);
	 BOOST_TEST(term14.key[0] == 1);
	 BOOST_TEST(term14.key[1] == 2);
	 BOOST_TEST(term14.key[2] == 3);
	 BOOST_TEST(term14.cf == -6.0);
 }

//void findTerm
 BOOST_AUTO_TEST_CASE(findTerm)
 {
	 Psym p1("p1");
	 Psym p2("p2");
	 Psym p3("p3");
	 VectorPsym polySym = { p1, p2, p3 };
	 boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	 SSeries::TermType::KeyType key01;
	 key01.resize(3);

	 // create a vector of 5 terms
	 std::vector<SSeries::TermType> sourceVector;
	 for (int i = 0; i < 5; ++i)
	 {
		 // we reuse the key
		 key01[0] = 1 * i;
		 key01[1] = 2 * i;
		 key01[2] = -2 * i;

		 SSeries::TermType::CfType value01(1.0 * i, polyOnlyArgs);
		 sourceVector.push_back(SSeries::TermType(value01, key01));
	 }

	 BOOST_TEST(sourceVector.size() == 5);  // just to verify we actually have 5 terms

	 SSeries bseries01; // the baseseries we populate
	 for (int i = 0; i < 5; ++i)
	 {
		 BOOST_CHECK_NO_THROW(bseries01.insert(sourceVector[i], polyOnlyArgs));
	 }

	 BOOST_TEST(bseries01.length() == 4); // all elements are different should have length of vector-1.
										  // index 0 constructs a null term(!) and is discarded	

	 // index 0 term shall not be in the series (a null term)
	 BOOST_TEST((bseries01.findTerm(sourceVector[0]) == bseries01.end()));

	 for (int i = 1; i < 5; ++i)
	 {
		 BOOST_TEST((*(bseries01.findTerm(sourceVector[i])) == sourceVector[i]));
	 }
 }


 //template <class Iterator, class ArgsTuple>
 //void insertRange(const Iterator&, const Iterator&, const ArgsTuple&);
 BOOST_AUTO_TEST_CASE(insertRange)
 {
	 Psym p1("p1");
	 Psym p2("p2");
	 Psym p3("p3");
	 VectorPsym polySym = { p1, p2, p3 };
	 boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	 SSeries::TermType::KeyType key01;
	 key01.resize(3);

	 // create a vector of 5 terms
	 std::vector<SSeries::TermType> sourceVector;
	 for (int i = 0; i < 5; ++i)
	 {
		 // we reuse the key
		 key01[0] = 1 * i;
		 key01[1] = 2 * i;
		 key01[2] = -2 * i;

		 SSeries::TermType::CfType value01(1.0*i, polyOnlyArgs);
		 sourceVector.push_back(SSeries::TermType(value01, key01));
	 }

	 BOOST_TEST(sourceVector.size() == 5);  // just to verify we actually have 5 terms

	 SSeries bseries01; // the baseseries we populate
	 BOOST_CHECK_NO_THROW(bseries01.insertRange(sourceVector.begin(), sourceVector.end(), polyOnlyArgs));
	 
	 BOOST_TEST(bseries01.length() == 4); // all elements are different should have length of vector-1.
										  // index 0 constructs a null term(!) and is discarded	

	 // index 0 term shall not be in the series (a null term)
	 BOOST_TEST((bseries01.findTerm(sourceVector[0]) == bseries01.end()));

	 for (int i = 1; i < 5; ++i)
	 {
		 BOOST_TEST((*(bseries01.findTerm(sourceVector[i])) == sourceVector[i]));
	 }
	 int teststop = 1;
 }

 //void baseSwap(Derived&);
 BOOST_AUTO_TEST_CASE(baseSwap)
 {
	 Psym p1("p1");
	 Psym p2("p2");
	 Psym p3("p3");
	 VectorPsym polySym = { p1, p2, p3 };
	 boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	 SSeries::TermType::KeyType key01;
	 SSeries::TermType::KeyType key02;
	 key01.resize(3);
	 key02.resize(3);
	 
	 // create a series of 4/5 Terms terms
	 SSeries bseries01; // the baseseries we populate
	 SSeries bseries02; // the baseseries we populate
	 for (int i = 0; i < 5; ++i)
	 {
		 // we reuse the key
		 key01[0] = 1 * i;
		 key01[1] = 2 * i;
		 key01[2] = -2 * i;

		 key02[0] = key01[0] + 1;
		 key02[1] = key01[1] + 1;
		 key02[2] = key01[2] + 1;

		 SSeries::TermType::CfType value01(1.0 * i, polyOnlyArgs);
		 SSeries::TermType::CfType value02(1.0 * i + 1, polyOnlyArgs);
		 bseries01.insert(SSeries::TermType(value01, key01), polyOnlyArgs);
		 bseries02.insert(SSeries::TermType(value02, key02), polyOnlyArgs);
	 }

	 BOOST_TEST(bseries01.length() == 4);
	 BOOST_TEST(bseries02.length() == 5);

	 bseries01.baseSwap(bseries02);
	 BOOST_TEST(bseries01.length() == 5);
	 BOOST_TEST(bseries02.length() == 4);
 }


 //size_type length()       const;
 BOOST_AUTO_TEST_CASE(length)
 {
     Psym p1("p1");
     VectorPsym polySym = { p1 };
     boost::tuple<VectorPsym> polyOnlyArgs(polySym);

     SSeries::TermType::KeyType key01;
     key01.resize(1);
     SSeries bseries01; // the baseseries we populate

     for (int i = 0; i < 5; ++i)
     {
         // we reuse the key
         key01[0] = 1 * i;

         SSeries::TermType::CfType value01(1.0 * i, polyOnlyArgs);
         bseries01.insert(SSeries::TermType(value01, key01), polyOnlyArgs);
     }

     BOOST_TEST(bseries01.length() == 4); // value i = 0 disappears!
 }


 //bool      empty()        const;
 BOOST_AUTO_TEST_CASE(emptySeries)
 {
     Psym p1("p1");
     VectorPsym polySym = { p1 };
     boost::tuple<VectorPsym> polyOnlyArgs(polySym);

     SSeries::TermType::KeyType key01;
     key01.resize(1);
     SSeries bseries01; // the baseseries we populate


     BOOST_TEST((bseries01.empty() == true)); // freshly constructed series

     SSeries bseries02;
     key01[0] = 1;
     SSeries::TermType::CfType value01(0.0, polyOnlyArgs);
     bseries01.insert(SSeries::TermType(value01, key01), polyOnlyArgs);

     BOOST_TEST(bseries02.empty() == true); // 0 added still keeps the series empty

     SSeries bseries03;
     key01[0] = 1;
     SSeries::TermType::CfType value02(1.0, polyOnlyArgs);
     bseries01.insert(SSeries::TermType(value02, key01), polyOnlyArgs);

     BOOST_TEST(bseries03.empty() == true); // 1 member added
 }


 //bool      isSingleCf()   const;
 BOOST_AUTO_TEST_CASE(isSingleCf)
 {
     Psym p1("p1");
     VectorPsym polySym = { p1 };
     boost::tuple<VectorPsym> polyOnlyArgs(polySym);

     {
         SSeries::TermType::KeyType key01;
         key01.resize(1);
         SSeries bseries01; // the baseseries we populate

         key01[0] = 0;
         SSeries::TermType::CfType value01(0.0, polyOnlyArgs);
         bseries01.insert(SSeries::TermType(value01, key01), polyOnlyArgs);
         BOOST_TEST(bseries01.isSingleCf() == false); // 
     }

     {
         SSeries::TermType::KeyType key02;
         key02.resize(1);
         SSeries bseries02;
         key02[0] = 0;
         SSeries::TermType::CfType value02(1.0, polyOnlyArgs);
         bseries02.insert(SSeries::TermType(value02, key02), polyOnlyArgs);
         BOOST_TEST(bseries02.isSingleCf() == true); // a constant
     }

     {
         SSeries::TermType::KeyType key03;
         key03.resize(1);
         SSeries bseries03;
         key03[0] = 1;
         SSeries::TermType::CfType value03(1.0, polyOnlyArgs);
         bseries03.insert(SSeries::TermType(value03, key03), polyOnlyArgs);
         BOOST_TEST(bseries03.isSingleCf() == false); // a constant
     }
 }

 //size_type atoms()        const;
 BOOST_AUTO_TEST_CASE(atoms) // not sure what it is good for?
 {
     Psym p1("p1");
     VectorPsym polySym = { p1 };
     boost::tuple<VectorPsym> polyOnlyArgs(polySym);

     {
         SSeries::TermType::KeyType key01;
         key01.resize(1);
         SSeries bseries01; // the baseseries we populate
         key01[0] = 0;
         SSeries::TermType::CfType value01(1.0, polyOnlyArgs);
         bseries01.insert(SSeries::TermType(value01, key01), polyOnlyArgs);
         BOOST_TEST(bseries01.atoms() == 1); //  ??
     }

     {
         SSeries::TermType::KeyType key02;
         key02.resize(1);
         SSeries bseries02; // the baseseries we populate
         key02[0] = 1;
         SSeries::TermType::CfType value02(1.0, polyOnlyArgs);
         bseries02.insert(SSeries::TermType(value02, key02), polyOnlyArgs);
         BOOST_TEST(bseries02.atoms() == 2); //  ??
     }

     {
         SSeries::TermType::KeyType key03;
         key03.resize(1);
         SSeries bseries03; // the baseseries we populate
         key03[0] = 1;
         SSeries::TermType::CfType value03(3.0, polyOnlyArgs);
         bseries03.insert(SSeries::TermType(value03, key03), polyOnlyArgs);
         BOOST_TEST(bseries03.atoms() == 0); //  ??
     }

 }


 // not sure what this is good for. requires series_multiplication and truncator
// I don't think it is used anywhere. except to define split in cf_series and named_series
// 
//template <class Series, class ArgsTuple>
//void baseSplit(std::vector<std::vector<Series> >&, const int n, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(baseSplit)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //
 // don't understand how that one really works.
 // test with Echelon level > 0
 //
 // create a vector of the series term. Basically use distributive law.
 // (x+2y+z)*cos(a) -> vector(x*cos, 2y*cos(a), z*cos(a))
 //
 //template <class ArgsTuple>
 //std::vector<TermType> flattenTerms(const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(flattenTerms)
 {
     BOOST_TEST_FAIL("No test implemented");
 }

 //
 // layouts apply to named series and describe the relation of the arguments to the key
 //
 //
 // applies a new layout to the terms of the BaseSeries. LayoutTuple is tuple (number of echonLevel members)
 // of vectors of pairs (bool, int). 
 // In BaseSeries we iterate over the the terms of the BaseSeries term container. 
 // for the cf as well as the key of the term we iterate over the size of the vector in the tuple and create a new
 // container of this size. The size of the current BasSeries term cf/key container has to be smaller or equal
 // to the size of the echelonLevel corresponding vector of the LayoutTuple. If the pair in the vector is 
 // (true,int) int gives an index into the current container and this index will be transfered to
 // new container at the current index position. Consequentially the int must be samller (it is an index) than the size of the current 
 // BaseSeries container. 
 //
 //template <class LayoutTuple, class ArgsTuple>
 //void applyLayoutToTerms(LayoutTuple const&, ArgsTuple const&, Derived&) const;
 BOOST_AUTO_TEST_CASE(applyLayoutToTerms)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //
 // test if we can trim some arguments from the argsTuple and terms. They can be trimmed if 
 // return the result as trimFlags(tuple of vector<bool>, false: trim i.e. the symbol in the key of the term is not used
 // The trimFlags vector has to be allocated and sized before calling trimTestTerms and has to agree with the sizes of the
 // key(s) used for the terms.
 // Used in cf_series and named_series.
 //
 //template <class TrimFlags>
 //void trimTestTerms(TrimFlags&) const;
 BOOST_AUTO_TEST_CASE(trimTestTerms)
 {
     BOOST_TEST_FAIL("No test implemented");
 }

 //
 // trim the keys according to the trimFlags from the terms and return the updated series in retval
 // The change of the argumentsTuple was done earlier before the use of trimTerms and the return series is setup
 // accordingly.
 //
 //template <class TrimFlags, class ArgsTuple>
 //void trimTerms(TrimFlags const&, ArgsTuple const&, Derived&) const;
 BOOST_AUTO_TEST_CASE(trimTerms)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //
 // To be clarified
 //
 //template <class RetSeries, class SubFunctor, class PosTuple, class SubCaches, class ArgsTuple>
 //RetSeries baseSub(const PosTuple&, SubCaches&, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(baseSub)
 {
     BOOST_TEST_FAIL("No test implemented");
 }

 //
 // first element of container
 //
 //const_iterator begin() const;
 BOOST_AUTO_TEST_CASE(beginSeries)
 {
     BOOST_TEST_FAIL("No test implemented");
 }

 //
 // end of container 
 //
 //const_iterator end() const;
 BOOST_AUTO_TEST_CASE(endSeries)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 /** @name Base comparison methods. */
 //template <class T>
 //bool baseEqualTo(const T&) const;
 BOOST_AUTO_TEST_CASE(baseEqualTo)
 {
     BOOST_TEST_FAIL("No test implemented");
 }

 //@}
 /** @name Base maths. */
 //@{
 //
 // calculate the norm of the series
 //
 //template <class ArgsTuple>
 //double baseNorm(const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(baseNorm)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //
 // evaluate series at value (a numeric value)
 // result (EvalType) can be double or complex<double> depending on series
 //
 //template <class ArgsTuple>
 //EvalType baseEval(const double, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(baseEval)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class T, class ArgsTuple>
 //Derived& baseAdd(const T&, const ArgsTuple&);
 // base addition of two series ??
 // can the ArgsTuple in the class and the explicit one in the method be different??
 // or has that to be checked outside of it???
 // what is ArgTuple actually good for? The symbols are nowhere used!!!!
 BOOST_AUTO_TEST_CASE(baseAdd)
 {
     Psym p1("p1");
     Psym p2("p2");
     VectorPsym polySym1 = { p1 };
     VectorPsym polySym2 = { p2 };
     VectorPsym polySym12 = { p1, p2 };
     boost::tuple<VectorPsym> polyOnlyArgs1(polySym1);
     boost::tuple<VectorPsym> polyOnlyArgs2(polySym2);
     boost::tuple<VectorPsym> polyOnlyArgs12(polySym12);

     // add two constants based on the same arguments
     {
         SSeries::TermType::KeyType key01;
         key01.resize(1);
         SSeries bseries01; // the baseseries we populate
         key01[0] = 0;
         SSeries::TermType::CfType value01(1.0, polyOnlyArgs1);
         bseries01.insert(SSeries::TermType(value01, key01), polyOnlyArgs1); // constant 1

         SSeries::TermType::KeyType key02;
         key02.resize(1);
         SSeries bseries02; // the baseseries we populate
         key02[0] = 0;
         SSeries::TermType::CfType value02(2.0, polyOnlyArgs1);
         bseries02.insert(SSeries::TermType(value02, key02), polyOnlyArgs1); // constant 2
         
         //add constants with same arguments p1
         bseries02.baseAdd(bseries01, polyOnlyArgs1);

         BOOST_TEST(bseries02.length() == 1);
         SSeries::TermType term02 = *(bseries02.begin());
         BOOST_TEST(term02.cf == 3.0);
         BOOST_TEST(term02.key[0] == 0);
         //BOOST_TEST(bseries02.)
     }

     // add two constants based on two different arguments
     {

     }

     // add a numerical constant to a constant series
     {

     }

     // add two non constant series with same key for same argument
     {

     }

     // add two non constant series with different keys for same argument
     {

     }

     // add two non constant series with different arguments
     {

     }
     
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class T, class ArgsTuple>
 //Derived& baseSubtract(const T&, const ArgsTuple&);
 BOOST_AUTO_TEST_CASE(baseSubtract)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //
 // requires series_multiplication to be usefull
 //
 //template <class T, class ArgsTuple>
 //Derived& baseMultBy(const T&, const ArgsTuple&);
 BOOST_AUTO_TEST_CASE(baseMultBy)
 {
     // in base_series_math.h
     // we have to distinguish multiplication with a coefficient and a multiplication with another 
     // series.
     // the series we multiply can have numerical coefficients only or they are series themselves!
     // selector is BaseSeriesMultiplySelector in base_series_mp.h
     //            
     // 1.     multiplication with constant coefficient(numerical)
     //         uses multiplyCoefficientsBy in in base_series_mp.h
     //         uses MultDivCoefficientsChecker in base_series_math.h
     // 1.1     multiplication with integer
     // 1.1.1     multiplication with 0
     // 1.1.2     mulitplication with unity 
     // 1.1.3     multiplication with other 
     //           uses multDivCoefficientsBy in base_series_math.h
     //             uses MultDivCoefficientsHelper in base_series_math.h
     //               uses multby in numerical_container.h or cf_series_math.h 
     // 1.2     multiplication with rational
     // 1.2.1     multiplication with 0
     // 1.2.2     mulitplication with unity 
     // 1.2.3     multiplication with other 
     //          uses multDivCoefficientsBy in base_series_math.h
     //             uses MultDivCoefficientsHelper in base_series_math.h

     // 1.3     multiplication with real
     // 1.3.1     multiplication with 0
     // 1.3.2     mulitplication with unity 
     // 1.3.3     multiplication with other 
     //           uses multDivCoefficientsBy in base_series_math.h
     //             uses MultDivCoefficientsHelper in base_series_math.h

     // 1.4     multiplication with complex
     // 1.4.1     multiplication with 0
     // 1.4.2     mulitplication with unity 
     // 1.4.3     multiplication with other 
     //           uses multDivCoefficientsBy in base_series_math.h 
     //             uses MultDivCoefficientsHelper in base_series_math.h



     // 2.  multiplication with series

     BOOST_TEST_FAIL("No test implemented");

 }


 //template <class T, class ArgsTuple>
 //Derived& baseDivideBy(const T&, const ArgsTuple&);
 BOOST_AUTO_TEST_CASE(baseDivideBy)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class ArgsTuple>
 //Derived baseInvert(const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(baseInvert)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class ArgsTuple>
 //Derived basePow(const double&, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(basePowDouble)
 {
     BOOST_TEST_FAIL("No test implemented");
 }

 //template <class ArgsTuple>
 //Derived basePow(const mp_rational&, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(basePowRational)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class ArgsTuple>
 //Derived naturalPower(const std::size_t, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(naturalPower)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class ArgsTuple>
 //Derived negativeIntegerPower(const int, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(negativeIntegerPower)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class ArgsTuple>
 //Derived realPower(const double, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(realPower)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class ArgsTuple>//
 //Derived rationalPower(const mp_rational&, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(rationalPower)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class ArgsTuple>
 //Derived baseRoot(const int, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(baseRoot)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class Series, class PosTuple, class ArgsTuple>
 //static void basePartial(const Derived&, Series&, const PosTuple&, const ArgsTuple&);
 BOOST_AUTO_TEST_CASE(basePartial)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class PosTuple, class ArgsTuple>
 //Derived basePartial(int, const PosTuple&, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(basePartial2)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class PosTuple, class ArgsTuple>
 //Derived basePartial(const PosTuple&, const ArgsTuple&) const;
 //@}
 BOOST_AUTO_TEST_CASE(basePartial3)
 {
     BOOST_TEST_FAIL("No test implemented");
 }

 /** @name Base output streaming methods. */
 //@{
 //template <class ArgsTuple>
 //void printTermsPlain(std::ostream&, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(printtermPlain)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class ArgsTuple>
 //void printTermsTEX(std::ostream&, const ArgsTuple&) const;
 BOOST_AUTO_TEST_CASE(printTermsTEX)
 {
     BOOST_TEST_FAIL("No test implemented");
 }


 //template <class ArgsTuple>
 //void printTermsPretty(std::ostream&, const ArgsTuple&) const;
 //@}
 BOOST_AUTO_TEST_CASE(printTermsPretty)
 {
     BOOST_TEST_FAIL("No test implemented");
 }



