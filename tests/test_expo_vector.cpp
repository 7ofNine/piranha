#define BOOST_TEST_MODULE Expo vector Test
#include "boost/test/included/unit_test.hpp"
#include "boost/test/tools/output_test_stream.hpp"

#include "piranha.h"


using namespace std;
using namespace piranha;
using namespace::manipulators;
using boost::test_tools::output_test_stream;
using namespace boost::unit_test;

namespace {

	using KeyType = ExpoVector<short int, 0>;// 0 is the position of the key i.e. echelon level
	using KeyRational = ExpoVector<mp_rational, 0>; // rational exponents are possible. when are they actually really used?
}

void setup() { PsymManager::clear(); }
BOOST_AUTO_TEST_CASE(simple_construction, *fixture(&setup))
{
	KeyType constructed;
	BOOST_TEST(constructed.size() == 0);
}

BOOST_AUTO_TEST_CASE(string_construction, *fixture(&setup))
{
	Psym p1("p1");
	Psym p2("p2");
	Psym p3("p3");
	Psym p4("p4");
	Psym p5("p5");
	VectorPsym polySym = { p1, p3, p5 };
	boost::tuple<VectorPsym> polyOnlyArgs(polySym);
	string values("1; 2;3"); // has deliberately a blank in it

	KeyType constructed3(values, polyOnlyArgs);
	BOOST_TEST(constructed3.size() == 3);
	BOOST_TEST((constructed3[0] == 1));
	BOOST_TEST((constructed3[1] == 2));
	BOOST_TEST((constructed3[2] == 3));


	// what does actually do
	KeyType constructed4(" 1", polyOnlyArgs);
	BOOST_TEST(constructed4.size() == 1);

	//BOOST_REQUIRE_THROW(( BaseExpoVector constructed5(" ", polyOnlyArgs), assertion_error)); ?
}

BOOST_AUTO_TEST_CASE(psym_construction, *fixture(&setup))
{
	// these are the same tests as for VectorKey
	Psym t1("t1");
	Psym t2("t2");

	VectorPsym polySym = { t1 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	{
		// Is that really used anywhere ?
		BOOST_CHECK_NO_THROW(KeyType tmp(t1, 0, polyOnlyArgs));
		KeyType tmp(t1, 0, polyOnlyArgs);
		BOOST_TEST(tmp.size() == 1);
		BOOST_TEST(tmp[0] == 1);
	}

	// symbol doesn't match
	BOOST_REQUIRE_THROW(KeyType tmp(t2, 0, polyOnlyArgs), assertion_error);

	{
		//position doesn't match, acts like an empty constructor
		KeyType tmp(t1, 2, polyOnlyArgs);
		BOOST_TEST(tmp.size() == 0);
	}

	// more than one symbol
	{
		VectorPsym polySym{ t1,t2 };
		boost::tuple<VectorPsym> polyOnlyArgs(polySym);
		BOOST_REQUIRE_THROW(KeyType tmp(t2, 0, polyOnlyArgs), assertion_error);
	}
}


BOOST_AUTO_TEST_CASE(multiply, *fixture(&setup))
{
	KeyType expo1;
	expo1.resize(3);
	expo1[0] = 1;
	expo1[1] = 2;
	expo1[2] = 3;

	KeyType expo2;
	expo2.resize(2);
	expo2[0] = 4;
	expo2[1] = 5;

	KeyType result;

	// multpily (i.e. add exponents)
	BOOST_CHECK_NO_THROW(expo1.multiply(expo2, result));
	BOOST_TEST(result.size() == 3);
	BOOST_TEST(result[0] == 5);
	BOOST_TEST(result[1] == 7);
	BOOST_TEST(result[2] == 3);

	// first factor must have more elements (why?)
	BOOST_CHECK_THROW(expo2.multiply(expo1, result), assertion_error);
}

BOOST_AUTO_TEST_CASE(printPlain, *fixture(&setup))
{
	// corresponds to printElements in vectorKey
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;

	// we have to have symbols
	Psym t1("t1");
	Psym t2("t2");
	Psym t3("t3");
	VectorPsym polySym{ t1, t2, t3 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	output_test_stream output;
	BOOST_CHECK_NO_THROW(temp.printPlain(output, polyOnlyArgs));
	BOOST_TEST(output.is_equal("1;2;3"));

	temp.invertSign();
	BOOST_CHECK_NO_THROW(temp.printPlain(output, polyOnlyArgs));
	BOOST_TEST(output.is_equal("-1;-2;-3"));
	
	// wrong number of args, they have to agree with the number of values in the vector
	VectorPsym polySym2{ t1,t3 };
	boost::tuple<VectorPsym> polyOnly2Args(polySym2);
	BOOST_CHECK_THROW(temp.printPlain(output, polyOnly2Args), assertion_error);

	VectorPsym polySym4{ t1,t3,t3,t3 };
	boost::tuple<VectorPsym> polyOnly4Args(polySym4);
	BOOST_CHECK_THROW(temp.printPlain(output, polyOnly4Args), assertion_error);
}

BOOST_AUTO_TEST_CASE(printPlainSorted, *fixture(&setup))
{
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	Psym t1("t1");
	Psym t2("t2");
	Psym t3("t3");
	VectorPsym polySym{ t1, t2, t3 };
	boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	output_test_stream output;
	std::vector<std::pair<bool, std::size_t>> sort = { {true, 2}, {true, 0}, {true, 1} };


	BOOST_CHECK_NO_THROW(temp.printPlainSorted(output, sort, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("     3      1      2 "));


	//sizes don't agree
	BOOST_CHECK_THROW(temp.printPlainSorted(output, std::vector<std::pair<bool, std::size_t>>(), polyOnlyArgs), assertion_error);

}

BOOST_AUTO_TEST_CASE(printPretty, *fixture(&setup))
{
	Psym t1("t1");
	Psym t2("t2");
	Psym t3("t3");
	VectorPsym polySym = { t1, t3, t2 };
	VectorPsym manySym = { t1,t1,t1,t1 };
	VectorPsym noSym;
	boost::tuple<VectorPsym> polyOnlyArgs(polySym);
	boost::tuple<VectorPsym> manyArgs(manySym);
	boost::tuple<VectorPsym> noArgs(noSym);

	KeyType term;
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
	BOOST_CHECK_NO_THROW(term.printPretty(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("t1*t3*t2"));

	// leading exponent is 0
	term[0] = 0;
	term[1] = 1;
	term[2] = 1;
	BOOST_CHECK_NO_THROW(term.printPretty(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("t3*t2"));

	// leading exponent is 0 and last is >1
	term[0] = 0;
	term[1] = 1;
	term[2] = 3;
	BOOST_CHECK_NO_THROW(term.printPretty(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("t3*t2**3"));

	// middle one is >1
	term[0] = 1;
	term[1] = 3;
	term[2] = 1;
	BOOST_CHECK_NO_THROW(term.printPretty(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("t1*t3**3*t2"));

	KeyRational tempr;
	tempr.resize(3);
	tempr[0] = 1;
	tempr[1] = mp_rational("2 / 3");
	tempr[2] = 0.5;
	BOOST_CHECK_NO_THROW(tempr.printPretty(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("t1*t3**(2/3)*t2**(1/2)"));

	tempr.resize(3);
	tempr[0] = 0;
	tempr[1] = mp_rational("2 / 3");
	tempr[2] = 0.5;
	BOOST_CHECK_NO_THROW(tempr.printPretty(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("t3**(2/3)*t2**(1/2)"));

	// more cases??
}

BOOST_AUTO_TEST_CASE(printTex, *fixture(&setup))
{
	Psym t1("t1");
	Psym t2("t2");
	Psym t3("t3");
	VectorPsym polySym = { t1, t3, t2 };
	VectorPsym manySym = { t1,t1,t1,t1 };
	VectorPsym noSym;
	boost::tuple<VectorPsym> polyOnlyArgs(polySym);
	boost::tuple<VectorPsym> manyArgs(manySym);
	boost::tuple<VectorPsym> noArgs(noSym);

	KeyType term;
	term.resize(3);
	term[0] = 1;
	term[1] = 1;
	term[2] = 1;
	output_test_stream output;

	// throw if sizes between ExpoVector and ArgsTuple component
	// are different
	BOOST_CHECK_THROW(term.printTex(output, noArgs), assertion_error);
	BOOST_CHECK_THROW(term.printTex(output, manyArgs), assertion_error);

	// everything is exponent 1
	BOOST_CHECK_NO_THROW(term.printTex(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal(" t1  t3  t2 ")); // two spaces between names

	// leading exponent is 0
	term[0] = 0;
	term[1] = 1;
	term[2] = 1;
	BOOST_CHECK_NO_THROW(term.printTex(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal(" t3  t2 "));

	// leading exponent is 0 and last is >1
	term[0] = 0;
	term[1] = 1;
	term[2] = 3;
	BOOST_CHECK_NO_THROW(term.printTex(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal(" t3  t2 ^{3}"));

	// middle one is >1
	term[0] = 1;
	term[1] = 3;
	term[2] = 1;
	BOOST_CHECK_NO_THROW(term.printTex(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal(" t1  t3 ^{3} t2 "));

	KeyRational tempr;
	tempr.resize(3);
	tempr[0] = 1;
	tempr[1] = mp_rational("2 / 3");
	tempr[2] = 0.5;
	BOOST_CHECK_NO_THROW(tempr.printTex(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal(" t1  t3 ^{\\frac{2}{3}} t2 ^{\\frac{1}{2}}")); // \\ doubled otherwise it acts like an escape character

	tempr.resize(3);
	tempr[0] = 0;
	tempr[1] = mp_rational("2 / 3");
	tempr[2] = 0.5;
	BOOST_CHECK_NO_THROW(tempr.printTex(output, polyOnlyArgs));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal(" t3 ^{\\frac{2}{3}} t2 ^{\\frac{1}{2}}"));

	// more cases?? neagtive numbers ??
}

BOOST_AUTO_TEST_CASE(evalTest, *fixture(&setup))
{
	Psym t1("t1", "-1.0;2.0");
	Psym t2("t2", "-3.0;4.0");
	Psym t3("t3", "5.0;6.0");
	VectorPsym polySym{ t1, t2, t3 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;

	BOOST_TEST(temp.eval(0.0, polyOnlyArgs) == -1125.0);
	BOOST_TEST(temp.eval(-1.0, polyOnlyArgs) == 147.0);
}

BOOST_AUTO_TEST_CASE(ignorable, *fixture(&setup))
{
	Psym t1("t1");
	Psym t2("t2");
	Psym t3("t3");
	VectorPsym polySym{ t1, t2, t3 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	
	// is never ignorable. Should we test that it actually has a size()>0
	// or does that directly correspond to all exponents are 0?
	BOOST_TEST(!temp.isIgnorable(polyOnlyArgs));
}

BOOST_AUTO_TEST_CASE(unity, *fixture(&setup))
{
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;

	// checks if all exponents are 0 i.e. the whole expo key == 1
	BOOST_TEST(!temp.isUnity());
	
	temp[0] = 0;
	temp[1] = 0;
	temp[2] = 0;
	BOOST_TEST(temp.isUnity());

	// what if key has no element?

}

BOOST_AUTO_TEST_CASE(comparison, *fixture(&setup))
{
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	KeyType temp1(temp);

	// they are equal i.e. fail
	BOOST_TEST(!(temp < temp1)); // (1,2,3) == (1,2,3)
	// all less than temp1 should pass
	temp[0]--;
	temp[1]--;
	temp[2]--;
	BOOST_TEST(temp < temp1);  // (0,1,2) < (1,2,3)
	// make one larger than temp1
	temp[1]++;
	BOOST_TEST(temp < temp1); //(0,2,2) < (1,2,3) 
	temp[0]++;
	BOOST_TEST(temp < temp1); //(1,2,2) < (1,2,3) 
	temp[1]++;
	BOOST_TEST(!(temp < temp1)); //(1,3,2) > (1,2,3) 
	// do we need more tests?

	// make unequal size
	temp.resize(5);
	BOOST_CHECK_THROW(temp.lexComparison(temp1), assertion_error);

}

BOOST_AUTO_TEST_CASE(norm_test, *fixture(&setup))
{
	Psym t1("t1","-1.0;2.0");
	Psym t2("t2","-3.0;4.0");
	Psym t3("t3","5.0;6.0");
	VectorPsym polySym{ t1, t2, t3 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;

	BOOST_TEST(temp.norm(polyOnlyArgs) == 1125);
}

BOOST_AUTO_TEST_CASE(hash_test, * fixture(&setup))
{
	//// not a real test. just to eecute the algorithm
	// this iidentical to wht is done for VectorKey
	KeyType temp;
	BOOST_TEST(temp.hash_value() == 0);

	temp.resize(3);
	temp[0] = -2;
	temp[1] = 3;
	temp[2] = 5;
	// Can we do something better here?
	BOOST_CHECK_NO_THROW(temp.hash_value());
}

BOOST_AUTO_TEST_CASE(degree, *fixture(&setup))
{
	Psym t1("t1", "-1.0;2.0", 3);
	Psym t2("t2", "-3.0;4.0", 2);
	Psym t3("t3", "5.0;6.0", 1);
	VectorPsym polySym{ t1, t2, t3 };
	
	KeyType temp;

	temp.resize(4);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	temp[3] = 4;

	BOOST_CHECK_THROW(temp.degree(polySym), assertion_error);

	temp = KeyType();
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;

	BOOST_TEST(temp.degree(polySym) == 10);
}


BOOST_AUTO_TEST_CASE(partialDegree, *fixture(&setup))
{
	Psym s1("n1", "1;2;3");
	Psym s2("n2");
	Psym s3("n3", "4.0;5.0;6.0", 7);
	Psym s4("n4");
	Psym s5("n5");

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

	VectorPsym polySym{ s1, s3, s5 };
	VectorPsym trigSym{ t1, t3, t5 };
	VectorPsym divSym{ d1, d3, d5 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);
	boost::tuple<VectorPsym, VectorPsym> polyTrigArgs(polySym, trigSym);
	boost::tuple<VectorPsym, VectorPsym, VectorPsym> polyTrigDivArgs(polySym, trigSym, divSym);
	boost::tuple<VectorPsym, VectorPsym, VectorPsym> allArgs{ {s1,s2,s3,s4,s5}, {t1,t2,t3,t4,t5}, {d1,d2,d3,d4,d5} };
	VectorPsym onePoly{ s3 };
	auto  posTuple = psyms2pos(onePoly, polyTrigDivArgs);

	KeyType temp;
	temp.resize(3);
	temp[0] = 3;
	temp[1] = 4;
	temp[2] = 5;

	// the consistency between argTuple and actual values is up to the user
	KeyType::DegreeType degree;
	BOOST_TEST((degree = temp.partialDegree(posTuple)) == 4);

	VectorPsym notExisting{ s4 };
	posTuple = psyms2pos(notExisting, polyTrigDivArgs);
	BOOST_TEST((degree = temp.partialDegree(posTuple)) == 0);

	VectorPsym threeplusAddes{ s5, s4,s1,s3 };
	posTuple = psyms2pos(threeplusAddes, polyTrigDivArgs);
	BOOST_TEST((degree = temp.partialDegree(posTuple)) == 12);
}

BOOST_AUTO_TEST_CASE(order, *fixture(&setup))
{
	// for ExpoVector this is defined tobe the same as degree 

	Psym t1("t1", "-1.0;2.0", 3);
	Psym t2("t2", "-3.0;4.0", 2);
	Psym t3("t3", "5.0;6.0", 1);
	VectorPsym polySym{ t1, t2, t3 };

	KeyType temp;

	temp.resize(4);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	temp[3] = 4;

	BOOST_CHECK_THROW(temp.order(polySym), assertion_error);

	temp = KeyType();
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;

	BOOST_TEST(temp.order(polySym) == 10);
}

BOOST_AUTO_TEST_CASE(partialOrder, *fixture(&setup))
{
	//for ExpoVector defined to be identical to partialDegree
	Psym s1("n1", "1;2;3");
	Psym s2("n2");
	Psym s3("n3", "4.0;5.0;6.0", 7);
	Psym s4("n4");
	Psym s5("n5");

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

	VectorPsym polySym{ s1, s3, s5 };
	VectorPsym trigSym{ t1, t3, t5 };
	VectorPsym divSym{ d1, d3, d5 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);
	boost::tuple<VectorPsym, VectorPsym> polyTrigArgs(polySym, trigSym);
	boost::tuple<VectorPsym, VectorPsym, VectorPsym> polyTrigDivArgs(polySym, trigSym, divSym);
	boost::tuple<VectorPsym, VectorPsym, VectorPsym> allArgs{ {s1,s2,s3,s4,s5}, {t1,t2,t3,t4,t5}, {d1,d2,d3,d4,d5} };
	VectorPsym onePoly{ s3 };
	auto  posTuple = psyms2pos(onePoly, polyTrigDivArgs);

	KeyType temp;
	temp.resize(3);
	temp[0] = 3;
	temp[1] = 4;
	temp[2] = 5;

	// the consistency between argTuple and actual values is up to the user
	KeyType::DegreeType degree;
	BOOST_TEST((degree = temp.partialOrder(posTuple)) == 4);

	VectorPsym notExisting{ s4 };
	posTuple = psyms2pos(notExisting, polyTrigDivArgs);
	BOOST_TEST((degree = temp.partialOrder(posTuple)) == 0);

	VectorPsym threeplusAddes{ s5, s4,s1,s3 };
	posTuple = psyms2pos(threeplusAddes, polyTrigDivArgs);
	BOOST_TEST((degree = temp.partialOrder(posTuple)) == 12);
}

BOOST_AUTO_TEST_CASE(partialDerivative, *fixture(&setup))
{
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

	Psym s1("n1","1;2;3");
	Psym s2("n2");
	Psym s3("n3", "4.0;5.0;6.0",7);
	Psym s4("n4");
	Psym s5("n5");

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

	VectorPsym polySym{ s1, s3, s5 };
	VectorPsym trigSym{ t1, t3, t5 };
	VectorPsym divSym { d1, d3, d5 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);
	boost::tuple<VectorPsym, VectorPsym> polyTrigArgs(polySym, trigSym);
	boost::tuple<VectorPsym, VectorPsym, VectorPsym> polyTrigDivArgs(polySym, trigSym, divSym);
	boost::tuple<VectorPsym, VectorPsym, VectorPsym> allArgs{ {s1,s2,s3,s4,s5}, {t1,t2,t3,t4,t5}, {d1,d2,d3,d4,d5} };
	VectorPsym onePoly{ s3 };
	//using PosTupleType = decltype(psyms2pos(onePoly, polyTrigDivArgs));
	auto  posTuple = psyms2pos(onePoly, polyTrigDivArgs);

	KeyType temp;
	temp.resize(3);
	temp[0] = 3;
	temp[1] = 4;
	temp[2] = 5;

	qpoly result;

	result = temp.partial<qpoly>(posTuple, polyTrigDivArgs);

	
	BOOST_TEST(result.length() == 1); // should only have a single element length is from BaseSeries.
	                                  // Most of the method should be from  BaseSeries.
	qpoly::TermType term;      // TermType is of class piranha::Monomial which is of type BaseTerm except for multiply method
	term = *(result.begin());  // get the only term
	KeyType resultKey = term.get<1>(); // get the key <1> gives the key
	auto coefficient = term.get<0>();  // everything else gives the coefficient only 0 and 1 are implemented!!
	BOOST_TEST(coefficient == 4);
	BOOST_TEST(resultKey[0] == 3);
	BOOST_TEST(resultKey[1] == 3);
	BOOST_TEST(resultKey[2] == 5);

	//check nothing has changed in the original symbol
	BOOST_TEST((s3.getTimeEval() == std::vector<double>{4.0, 5.0, 6.0}));
	BOOST_TEST(s3.order() = 7);

	VectorPsym twoSymbols{ s3, s5 };
	// only one symbol is allowed (partial derivative for one variable only)
	auto posTuple2 = psyms2pos(twoSymbols, polyTrigDivArgs);
	BOOST_CHECK_THROW(temp.partial<qpoly>(posTuple2, polyTrigDivArgs), assertion_error);


	// symbol outside the range of the vectorKey
	// the index is taken from where it is int ArgType
	VectorPsym outsideSymbol{ s5 };
	auto posTupleOut = psyms2pos(outsideSymbol, allArgs);
	BOOST_CHECK_THROW(temp.partial<qpoly>(posTupleOut, polyTrigDivArgs), assertion_error);

	// symbol not present
	//VectorPsym notPresent{ s4 };
	//auto noSymbol = psyms2pos(notPresent, polyTrigDivArgs);
	////nothing should happen. derivative for a symbol that doesn't exist. The result is 0;
	// this is currently not the case. The result has 0 length which mathematically is incorrect
	// TODO: FIX IT
	//result = temp.partial<qpoly>(noSymbol, polyTrigDivArgs);
	//BOOST_TEST(result.length() == 1);


}

BOOST_AUTO_TEST_CASE(power, *fixture(&setup))
{
	Psym s1("n1", "1;2;3");
	Psym s2("n2");
	Psym s3("n3", "4.0;5.0;6.0", 7);
	VectorPsym polySym{ s1, s3};
	boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	KeyType temp;
	temp.resize(3);
	temp[0] = 0;
	temp[1] = 2;
	temp[2] = 3;
	
	// integer power
	KeyType newTemp = temp.pow(2, polyOnlyArgs);
	BOOST_TEST(newTemp[0] == 0);
	BOOST_TEST(newTemp[1] == 4);
	BOOST_TEST(newTemp[2] == 6);


	// rational power on integer exponents should throw
	mp_rational power(1, 3);
	BOOST_CHECK_THROW(temp.pow(power, polyOnlyArgs), value_error);

	using KeyRational = ExpoVector<mp_rational, 0>;
	KeyRational tempr;
	tempr.resize(1);
	tempr[0] = mp_rational(2, 3);
	KeyRational resultr;
	resultr = tempr.pow(power, polyOnlyArgs);
	BOOST_TEST(resultr[0] == mp_rational("2/9"));

	// double on integer should throw
	BOOST_CHECK_THROW(temp.pow(1.5, polyOnlyArgs), value_error);

	// double on rational transformed to rational on rational
	resultr = tempr.pow(1.5, polyOnlyArgs);
	BOOST_TEST(resultr[0] == mp_rational(1, 1));

	// double only on unity (i.e. does nothing)
	temp[0] = 1;
	temp[1] = 1;
	temp[2] = 1;

	BOOST_CHECK_THROW(temp.pow(1.5, polyOnlyArgs),value_error);
	temp[0] = 0;
	temp[1] = 0;
	temp[2] = 0;
	newTemp = temp.pow(1.5, polyOnlyArgs);
	BOOST_TEST(newTemp[0] == 0);
	BOOST_TEST(newTemp[1] == 0);
	BOOST_TEST(newTemp[2] == 0);
}

BOOST_AUTO_TEST_CASE(A_substitute, *fixture(&setup))
{
	// how to implement a test? requires a series
	//should we wait until we understand where it is used and how to actually test it
	// possibly only used trhough common_poisson_series_toolbox
	// and named_series defines substitute
	//There are many layers involved before we get here. No understanding,yet.
	BOOST_FAIL("Needs investigation to implement");
}

BOOST_AUTO_TEST_CASE(A_eiSubstitute, *fixture(&setup))
{
	// is that actually anywhere used and functional??
	//should we wait until we understand where it is used and how to actually test it
	// possibly only used trhough common_poisson_series_toolbox
	BOOST_FAIL("Needs investigation to implement");
}

