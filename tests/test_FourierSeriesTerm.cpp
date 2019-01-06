#define BOOST_TEST_MODULE BaseSeries Test
#include "boost/test/included/unit_test.hpp"

#include "piranha.h"

// This also tests BaseTerm because of CRTP

using namespace std;
using namespace piranha;

// A FourierSeriesTerm has 4 template parameters 
// Cf: for coefficients which should be a piranha *cf class
// Trig: a trigonometric key
// Separator: for writing and reading purposes (typical '|')
// Allocator: typical std::allocator but we can use diffrent ones e.g. CountingAllocator
namespace {
	using SimpleFourierTerm   = FourierSeriesTerm<double, TrigVector<short, 0>, '|', std::allocator<char>>;// I can create the type but would that work??
	using SimpleCfFourierTerm = FourierSeriesTerm<double_cf, TrigVector<short, 0>, '|', std::allocator<char>>;
}


BOOST_AUTO_TEST_CASE(construction_test)
{
	Psym t1("t1");
	Psym t2("t2");
	Psym t3("t3");

	VectorPsym trigSymbols = { t1, t3 };
	boost::tuple<VectorPsym> trigOnlyArgs(trigSymbols);

	SimpleFourierTerm fourierTerm;
	//What to check?


	// why is argtuples of interest
	TrigVector<short, 0> testVector("1; -2;c", trigOnlyArgs);
	SimpleFourierTerm fourierTermfromCfKey(1.1, testVector);

	SimpleCfFourierTerm fourierCfTerm(double_cf(-2.2, trigOnlyArgs), testVector);

	//THis does not compile!
	// no coefficient given
	//SimpleFourierTerm onlyArguments("1.1|-1; 2;s", trigOnlyArgs);
	//SimpleFourierTerm onlyArgumentsCf("2.2|-1; 2;s", trigOnlyArgs);
	//coefficent given
	// 
	//SimpleFourierTerm withCfandArguments("1.1|-1; 2;s", trigOnlyArgs);
	//SimpleFourierTerm withCfandArgumentsCf("2.2|-1; 2;s", trigOnlyArgs);



	//BOOST_TEST_FAIL("No test implemented");

}

BOOST_AUTO_TEST_CASE(test_multiply)
{
	// doesn't compile, yet where is the problem??
	//Psym t1("t1");
	//Psym t2("t2");
	//Psym t3("t3");
	//Psym t4("t4");
	//Psym t5("t5");

	//VectorPsym trigSymbolsFactor1 = { t1, t3 };
	//boost::tuple<VectorPsym> argsFactor1(trigSymbolsFactor1);

	//VectorPsym trigSymbolsFactor2 = { t2, t3 };
	//boost::tuple<VectorPsym> argsFactor2(trigSymbolsFactor1);

	//TrigVector<short, 0> vectorFactorCos1("1; -2;c", trigSymbolsFactor1);
	//TrigVector<short, 0> vectorFactorSin2("-3; 5;s", trigSymbolsFactor1);

	//TrigVector<short, 0> vectorFactorCos3("7; -11;c",  trigSymbolsFactor1);
	//TrigVector<short, 0> vectorFactorSin4("-13; 17;s", trigSymbolsFactor1);

	//// there are four possible results
	//SimpleCfFourierTerm::multiplication_result result1;
	//SimpleCfFourierTerm::multiplication_result result2;
	//SimpleCfFourierTerm::multiplication_result result3;
	//SimpleCfFourierTerm::multiplication_result result4;

	//SimpleCfFourierTerm factor1(double_cf(1.0, trigSymbolsFactor1), vectorFactorCos1);
	//SimpleCfFourierTerm factor2(double_cf(1.0, trigSymbolsFactor1), vectorFactorCos3);
	//SimpleCfFourierTerm factor3(double_cf(1.0, trigSymbolsFactor1), vectorFactorSin2);
	//SimpleCfFourierTerm factor4(double_cf(1.0, trigSymbolsFactor1), vectorFactorSin4);

	//// cos*cos
	//SimpleCfFourierTerm::multiply(factor1, factor2, result1, argsFactor1);
	//// sin * sin
	//SimpleCfFourierTerm::multiply(factor3, factor4, result2, argsFactor1);
	//// sin * cos
	//SimpleCfFourierTerm::multiply(factor4, factor1, result3, argsFactor1);
	////cos * sin
	//SimpleCfFourierTerm::multiply(factor2, factor3, result3, argsFactor1);

}
