#define BOOST_TEST_MODULE Psym Test
#include "boost/test/included/unit_test.hpp"
#include "boost/test/output_test_stream.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

#include <string>
#include <set>

#include "piranha.h"

using namespace std;
using namespace piranha;
using boost::test_tools::output_test_stream;

namespace {
	static int const DEFAULT_ORDER = 1;
	static const string name1 = "name1";
	static const string name2 = "name2";
	static const string name3 = "name3";
	static const string name4 = "name4";
	static const string name5 = "name5";
	static const string name6 = "name6";
	static const string name7 = "name7";
	static const string name8 = "name8";
	static const string name9 = "name9";
}

BOOST_AUTO_TEST_CASE(construction_test)
{
	
	// construct with name only
	Psym symbol1(name1);
	BOOST_TEST(symbol1.getName() == "name1");
	BOOST_TEST(symbol1.getTimeEval().empty());
	BOOST_TEST(symbol1.order() == DEFAULT_ORDER);

	// construct with name and order
	Psym symbol2(name2, 3);
	BOOST_TEST(symbol2.getName() == "name2");
	BOOST_TEST(symbol2.getTimeEval().empty());
	BOOST_TEST(symbol2.order() == 3);

	// construct with name and time evaluation vector
	vector<double> timeEval{ 1.0, -1.1, +2.2, 3.3 };
	// the string deliberately has blanks in the middle, to check string normalization before casting
	string const timeEvalString("1.0; -1.1;+2.2 ; 3.3 "); // can we get to the separator ";" programatically ?
	Psym symbol3(name3, timeEval);
	BOOST_TEST(symbol3.getName() == name3);
	BOOST_TEST(symbol3.getTimeEval() == timeEval);
	BOOST_TEST(symbol3.order() = DEFAULT_ORDER);

	// construct with name and time evaluation vector and order
	unsigned int const testOrder = 99;
	Psym symbol4(name4, timeEval, testOrder);
	BOOST_TEST(symbol4.getName() == name4);
	BOOST_TEST(symbol4.getTimeEval() == timeEval);
	BOOST_TEST(symbol4.order() == testOrder);

	//construct with name and time evaluation string
	Psym symbol5(name5, timeEvalString);
	BOOST_TEST(symbol5.getName() == name5);
	BOOST_TEST(symbol5.getTimeEval() == timeEval);
	BOOST_TEST(symbol5.order() == DEFAULT_ORDER);


	////construct with name and time evaluation string and order
	Psym symbol6(name6, timeEvalString, testOrder);
	BOOST_TEST(symbol6.getName() == name6);
	BOOST_TEST(symbol6.getTimeEval() == timeEval);
	BOOST_TEST(symbol6.order() == testOrder);

	//construct with name and empty time evaluation string and order
	string const emptyTimeEvalString("");
	Psym symbol7(name7, emptyTimeEvalString, testOrder);
	BOOST_TEST(symbol7.getName() == name7);
	BOOST_TEST(symbol7.getTimeEval().empty());
	BOOST_TEST(symbol7.order() == testOrder);

	// construct the same name twice
	Psym symbol8(name8);
	size_t const oldSize = Psym::list().size();// preserve list size

	BOOST_TEST(symbol8.getName() == name8);
	BOOST_TEST(symbol8.getTimeEval().empty());
	BOOST_TEST(symbol8.order() == DEFAULT_ORDER);
	// should have the updated time evaluation polynomial and order
	Psym symbolx(name8, timeEvalString, testOrder);
	BOOST_TEST(symbolx.getName() == name8);
	BOOST_TEST(symbolx.getTimeEval() == timeEval);
	BOOST_TEST(symbolx.order() == testOrder);
	// and symbolx should be the same as symbol8
	BOOST_TEST(symbolx.getName() == symbol8.getName());
	BOOST_TEST(symbolx.getTimeEval() == symbol8.getTimeEval());
	BOOST_TEST(symbolx.order() == symbol8.order());
	BOOST_TEST((symbolx == symbol8));
	// symbol list size should not have changed
	BOOST_TEST(Psym::list().size() == oldSize);

	// create with same name (other parameters are kept)
	// the original behaviour
	Psym symbolu(name8);
	BOOST_TEST((symbolu == symbol8));
	BOOST_TEST((symbolu.getTimeEval() == symbol8.getTimeEval()));
	BOOST_TEST((symbolu.order() == symbol8.order()));


	// do double with name and order
	Psym symbolz(name8, 3);
	BOOST_TEST(symbolz.getName() == name8);
	BOOST_TEST((symbolz.getTimeEval() == symbol8.getTimeEval()));
	BOOST_TEST(symbolz.order() == 3);
	BOOST_TEST(Psym::list().size() == oldSize);

	// construct with name and single value (i.e.  a constant value of order 0 in the time evaluation)
	Psym symbol9(name9, 1.5);
	BOOST_TEST(symbol9.getName() == name9);
	BOOST_TEST(symbol9.getTimeEval().size() == 1);
	BOOST_TEST(symbol9.getTimeEval()[0] == 1.5);
	BOOST_TEST(symbol9.order() == DEFAULT_ORDER);
}

BOOST_AUTO_TEST_CASE(operator_test)
{
	// construct with name1
	Psym symbol1(name1);
	BOOST_TEST(symbol1.getName() == name1);

	// construct with name2 
	Psym symbol2(name2, 3);
	BOOST_TEST(symbol2.getName() == name2);

	//comparison
	BOOST_TEST(!(symbol2 < symbol1));

	BOOST_TEST((symbol1 < symbol2));

	BOOST_TEST((symbol1 != symbol2));

	Psym symbol(name1);
	BOOST_TEST((symbol1 == symbol));
}

BOOST_AUTO_TEST_CASE(setter_getter_test)
{
	//int const defaultOrder = 1; // the default order set in constructor
	Psym symbol1(name1);
	int const testOrder = 99;

	BOOST_TEST(symbol1.order() == DEFAULT_ORDER);

	// set new order
	symbol1.setOrder(testOrder);
	BOOST_TEST(symbol1.order() = testOrder);

	// test time evaluation polynomial
	// set new time evaluation vector via constructor
	vector<double> const testTimeEval = { 0.0, 1.0, 2.0 };

	Psym symbol2(name2, testTimeEval);
	BOOST_TEST((symbol2.getTimeEval() == testTimeEval));
	//evaluate time dependent value
	BOOST_TEST(symbol2.eval(0.0) == 0.0);
	BOOST_TEST(symbol2.eval(1.0) == 3.0);
	BOOST_TEST(symbol2.eval(-1.0) == 1.0);


	// set new time evaluation vector explicitly
	Psym symbol3(name3);
	//first it is empty
	BOOST_TEST((symbol3.getTimeEval().empty()));
	symbol3.setTimeEval(testTimeEval);
	//now it is set
	BOOST_TEST((symbol3.getTimeEval() == testTimeEval));
	//evaluate time dependent value
	BOOST_TEST(symbol3.eval(0.0) == 0.0);
	BOOST_TEST(symbol3.eval(1.0) == 3.0);
	BOOST_TEST(symbol3.eval(-1.0) == 1.0);
}

BOOST_AUTO_TEST_CASE(list_test)
{
	// at the beginning the list is empty
	BOOST_TEST(Psym::list().empty());

	// add 5 names to the list the primitive way
	Psym s1("n1");
	Psym s2("n2");
	Psym s3("n3");
	Psym s4("n4");
	Psym s5("n5");
	BOOST_TEST((Psym::list().size() == 5));
	
	//check list
	// transform vector into set. THe list is actually a vector
	auto const symbolList = Psym::list();
	set<Psym> testSet(begin(symbolList), end(symbolList)); // somehow it crashes when I do it directly from the Psym.list()
	// nothing should be lost or added, The Psym list is supposed to contain
	// unique items
	BOOST_TEST((Psym::list().size() == testSet.size()));
	// each element should now be in there once
	BOOST_TEST((testSet.count(s1) == 1));
	BOOST_TEST((testSet.count(s2) == 1));
	BOOST_TEST((testSet.count(s3) == 1));
	BOOST_TEST((testSet.count(s4) == 1));
	BOOST_TEST((testSet.count(s5) == 1));
	// now all elements that we expected are in the list and no additional ones

	// test for names2psyms
	vector<string> testVector1 = { "n4", "n1", "n5" };
	VectorPsym resultVector = names2psyms(testVector1); // transform names to a corresponding vector of psyms
	BOOST_TEST((resultVector[0] == s4));
	BOOST_TEST((resultVector[1] == s1));
	BOOST_TEST((resultVector[2] == s5));
	BOOST_TEST((resultVector.size() == 3));

	// test for psyms2pos
	// Add some pseudo trigonometric and divisor names (level 3 is not practically used right now)
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

	VectorPsym polySym = { s1, s3, s5 };
	VectorPsym trigSym = { t1, t3, t5 };
	VectorPsym divSym = { d1, d3, d5 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);
	boost::tuple<VectorPsym, VectorPsym> polyTrigArgs(polySym, trigSym);
	boost::tuple<VectorPsym, VectorPsym, VectorPsym> polyTrigDivArgs(polySym, trigSym, divSym);

	// there seems to be no need for these Ntuples to be constructed explicitly with symbols pre-set.
	// they are instantiation created an arguments explixitly set. we'll have to test that with named series
	// we construct them here for trial purposes 
	//NTuple<VectorPsym, 1> tpolyOnlyArgs;
	// tpolyOnlyArgs.get<0>().push_back(s1); ?? why would get<N> not work here while it works in piranha code?
	//piranha::NTuple<piranha::VectorPsym, 2> tpolytrigArgs;
	//piranha::NTuple<piranha::VectorPsym, 3> tpolytrigdivArgs;

	// the return type is a mouthfull, using auto
	// NTuple<std::vector<std::pair<bool, std::size_t> >, boost::tuples::length<ArgsTuple>::value>::Type
	// get<N> with incorrect values does not compile(!) they are templated.
	VectorPsym onePoly = { s3 };
	auto onePolyPolyResult = psyms2pos(onePoly, polyOnlyArgs);
	BOOST_TEST((onePolyPolyResult.get<0>().size() == 1));
	BOOST_TEST((onePolyPolyResult.get<0>()[0].first == true));
	BOOST_TEST((onePolyPolyResult.get<0>()[0].second == 1));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing

	auto onePolyDivResult  = psyms2pos(onePoly, polyTrigDivArgs);
	BOOST_TEST((onePolyDivResult.get<0>().size() == 1));
	BOOST_TEST((onePolyDivResult.get<0>()[0].first == true));
	BOOST_TEST((onePolyDivResult.get<0>()[0].second == 1));
	BOOST_TEST((onePolyDivResult.get<1>().size() == 1));
	BOOST_TEST((onePolyDivResult.get<1>()[0].first == false));
	BOOST_TEST((onePolyDivResult.get<1>()[0].second == 0));
	BOOST_TEST((onePolyDivResult.get<2>().size() == 1));
	BOOST_TEST((onePolyDivResult.get<2>()[0].first == false));
	BOOST_TEST((onePolyDivResult.get<2>()[0].second == 0));

	VectorPsym onePolynotPresent = { s2 };
	auto onePolyNotPresentResult = psyms2pos(onePolynotPresent, polyOnlyArgs);
	BOOST_TEST((onePolyNotPresentResult.get<0>().size() == 1));
	BOOST_TEST((onePolyNotPresentResult.get<0>()[0].first == false));
	BOOST_TEST((onePolyNotPresentResult.get<0>()[0].second == 0));

	auto onePolyDivNotPresentResult = psyms2pos(onePolynotPresent, polyTrigDivArgs);
	BOOST_TEST((onePolyDivNotPresentResult.get<0>().size() == 1));
	BOOST_TEST((onePolyDivNotPresentResult.get<0>()[0].first == false));
	BOOST_TEST((onePolyDivNotPresentResult.get<0>()[0].second == 0));
	BOOST_TEST((onePolyDivNotPresentResult.get<1>().size() == 1));
	BOOST_TEST((onePolyDivNotPresentResult.get<1>()[0].first == false));
	BOOST_TEST((onePolyDivNotPresentResult.get<1>()[0].second == 0));
	BOOST_TEST((onePolyDivNotPresentResult.get<2>().size() == 1));
	BOOST_TEST((onePolyDivNotPresentResult.get<2>()[0].first == false));
	BOOST_TEST((onePolyDivNotPresentResult.get<2>()[0].second == 0));

	// test double entry exception
	VectorPsym doubleEntry = { s3, s3 };
	// this should blow up
	BOOST_REQUIRE_THROW(psyms2pos(doubleEntry, polyOnlyArgs), assertion_error);

	// three polynomial entries in order
	VectorPsym triplePoly = { s1, s3, s5 };
	auto triplePolyPolyResult = psyms2pos(triplePoly, polyOnlyArgs);
	BOOST_TEST((triplePolyPolyResult.get<0>().size() == 3));
	BOOST_TEST((triplePolyPolyResult.get<0>()[0].first == true));
	BOOST_TEST((triplePolyPolyResult.get<0>()[0].second == 0));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing
	BOOST_TEST((triplePolyPolyResult.get<0>()[1].first == true));
	BOOST_TEST((triplePolyPolyResult.get<0>()[1].second == 1));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing
	BOOST_TEST((triplePolyPolyResult.get<0>()[2].first == true));
	BOOST_TEST((triplePolyPolyResult.get<0>()[2].second == 2));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing

		// three polynomial entries out of order
	VectorPsym triplePolyOrder = { s3, s5, s1 };
	auto triplePolyPolyOrderResult = psyms2pos(triplePolyOrder, polyOnlyArgs);
	BOOST_TEST((triplePolyPolyOrderResult.get<0>().size() == 3));
	BOOST_TEST((triplePolyPolyOrderResult.get<0>()[0].first == true));
	BOOST_TEST((triplePolyPolyOrderResult.get<0>()[0].second == 1));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing
	BOOST_TEST((triplePolyPolyOrderResult.get<0>()[1].first == true));
	BOOST_TEST((triplePolyPolyOrderResult.get<0>()[1].second == 2));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing
	BOOST_TEST((triplePolyPolyOrderResult.get<0>()[2].first == true));
	BOOST_TEST((triplePolyPolyOrderResult.get<0>()[2].second == 0));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing

	// mixed one of everything
	VectorPsym mixedSingleOrder = { s3, t5, d1 };
	auto mixedSingleOrderResult = psyms2pos(mixedSingleOrder, polyTrigDivArgs);
	BOOST_TEST((mixedSingleOrderResult.get<0>().size() == 3));
	BOOST_TEST((mixedSingleOrderResult.get<0>()[0].first == true));
	BOOST_TEST((mixedSingleOrderResult.get<0>()[0].second == 1));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing
	BOOST_TEST((mixedSingleOrderResult.get<0>()[1].first == false));
	BOOST_TEST((mixedSingleOrderResult.get<0>()[1].second == 0));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing
	BOOST_TEST((mixedSingleOrderResult.get<0>()[2].first == false));
	BOOST_TEST((mixedSingleOrderResult.get<0>()[2].second == 0));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing

	//trig arguments
	BOOST_TEST((mixedSingleOrderResult.get<1>().size() == 3));
	BOOST_TEST((mixedSingleOrderResult.get<1>()[0].first == false));
	BOOST_TEST((mixedSingleOrderResult.get<1>()[0].second == 0));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing
	BOOST_TEST((mixedSingleOrderResult.get<1>()[1].first == true));
	BOOST_TEST((mixedSingleOrderResult.get<1>()[1].second == 2));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing
	BOOST_TEST((mixedSingleOrderResult.get<1>()[2].first == false));
	BOOST_TEST((mixedSingleOrderResult.get<1>()[2].second == 0));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing

	// div arguments
	BOOST_TEST((mixedSingleOrderResult.get<2>().size() == 3));
	BOOST_TEST((mixedSingleOrderResult.get<2>()[0].first == false));
	BOOST_TEST((mixedSingleOrderResult.get<2>()[0].second == 0));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing
	BOOST_TEST((mixedSingleOrderResult.get<2>()[1].first == false));
	BOOST_TEST((mixedSingleOrderResult.get<2>()[1].second == 0));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing
	BOOST_TEST((mixedSingleOrderResult.get<2>()[2].first == true));
	BOOST_TEST((mixedSingleOrderResult.get<2>()[2].second == 0));  // the index into Psymanager. Sensitive to changes. Not really a good way of programing

	// mixed complete         0   1   2   3   4   5   6   7   8   (index in result vectors)
	VectorPsym mixedOrder = { s3, t5, s2, t1, d1, d5, t2, s1, d3};
	auto mixedOrderResult = psyms2pos(mixedOrder, polyTrigDivArgs);
	BOOST_TEST((mixedOrderResult.get<0>().size() == 9));
	BOOST_TEST((mixedOrderResult.get<0>()[0].first == true));
	BOOST_TEST((mixedOrderResult.get<0>()[0].second == 1));
	BOOST_TEST((mixedOrderResult.get<0>()[1].first == false));
	BOOST_TEST((mixedOrderResult.get<0>()[1].second == 0));
	BOOST_TEST((mixedOrderResult.get<0>()[2].first == false));
	BOOST_TEST((mixedOrderResult.get<0>()[2].second == 0));
	BOOST_TEST((mixedOrderResult.get<0>()[3].first == false));
	BOOST_TEST((mixedOrderResult.get<0>()[3].second == 0));
	BOOST_TEST((mixedOrderResult.get<0>()[4].first == false));
	BOOST_TEST((mixedOrderResult.get<0>()[4].second == 0));
	BOOST_TEST((mixedOrderResult.get<0>()[5].first == false));
	BOOST_TEST((mixedOrderResult.get<0>()[5].second == 0));
	BOOST_TEST((mixedOrderResult.get<0>()[6].first == false));
	BOOST_TEST((mixedOrderResult.get<0>()[6].second == 0));
	BOOST_TEST((mixedOrderResult.get<0>()[7].first == true));
	BOOST_TEST((mixedOrderResult.get<0>()[7].second == 0));
	BOOST_TEST((mixedOrderResult.get<0>()[8].first == false));
	BOOST_TEST((mixedOrderResult.get<0>()[8].second == 0));

	// trig arguments
	BOOST_TEST((mixedOrderResult.get<1>().size() == 9));
	BOOST_TEST((mixedOrderResult.get<1>()[0].first == false));
	BOOST_TEST((mixedOrderResult.get<1>()[0].second == 0));
	BOOST_TEST((mixedOrderResult.get<1>()[1].first == true));
	BOOST_TEST((mixedOrderResult.get<1>()[1].second == 2));
	BOOST_TEST((mixedOrderResult.get<1>()[2].first == false));
	BOOST_TEST((mixedOrderResult.get<1>()[2].second == 0));
	BOOST_TEST((mixedOrderResult.get<1>()[3].first == true));
	BOOST_TEST((mixedOrderResult.get<1>()[3].second == 0));
	BOOST_TEST((mixedOrderResult.get<1>()[4].first == false));
	BOOST_TEST((mixedOrderResult.get<1>()[4].second == 0));
	BOOST_TEST((mixedOrderResult.get<1>()[5].first == false));
	BOOST_TEST((mixedOrderResult.get<1>()[5].second == 0));
	BOOST_TEST((mixedOrderResult.get<1>()[6].first == false));
	BOOST_TEST((mixedOrderResult.get<1>()[6].second == 0));
	BOOST_TEST((mixedOrderResult.get<1>()[7].first == false));
	BOOST_TEST((mixedOrderResult.get<1>()[7].second == 0));
	BOOST_TEST((mixedOrderResult.get<1>()[8].first == false));
	BOOST_TEST((mixedOrderResult.get<1>()[8].second == 0));

	// div arguments
	// mixed complete            0   1   2   3   4   5   6   7   8   (index in result vectors)
	// VectorPsym mixedOrder = { s3, t5, s2, t1, d1, d5, t2, s1, d3 };
	BOOST_TEST((mixedOrderResult.get<2>().size() == 9));
	BOOST_TEST((mixedOrderResult.get<2>()[0].first == false));
	BOOST_TEST((mixedOrderResult.get<2>()[0].second == 0));
	BOOST_TEST((mixedOrderResult.get<2>()[1].first == false));
	BOOST_TEST((mixedOrderResult.get<2>()[1].second == 0));
	BOOST_TEST((mixedOrderResult.get<2>()[2].first == false));
	BOOST_TEST((mixedOrderResult.get<2>()[2].second == 0));
	BOOST_TEST((mixedOrderResult.get<2>()[3].first == false));
	BOOST_TEST((mixedOrderResult.get<2>()[3].second == 0));
	BOOST_TEST((mixedOrderResult.get<2>()[4].first == true));
	BOOST_TEST((mixedOrderResult.get<2>()[4].second == 0));
	BOOST_TEST((mixedOrderResult.get<2>()[5].first == true));
	BOOST_TEST((mixedOrderResult.get<2>()[5].second == 2));
	BOOST_TEST((mixedOrderResult.get<2>()[6].first == false));
	BOOST_TEST((mixedOrderResult.get<2>()[6].second == 0));
	BOOST_TEST((mixedOrderResult.get<2>()[7].first == false));
	BOOST_TEST((mixedOrderResult.get<2>()[7].second == 0));
	BOOST_TEST((mixedOrderResult.get<2>()[8].first == true));
	BOOST_TEST((mixedOrderResult.get<2>()[8].second == 1));



}

BOOST_AUTO_TEST_CASE(output_test)
{
	output_test_stream output;

	Psym symbol1(name1);
//	output_test_stream().swap(output); // clean output
	symbol1.print(output);
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("name=name1\ntime_eval=\norder=1\n"));

	// construct with name and order
	Psym symbol2(name2, 3);
	symbol2.print(output);
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(symbol2.getTimeEval().empty());
	BOOST_TEST(output.is_equal("name=name2\ntime_eval=\norder=3\n"));


	// construct with name and time evaluation vector
	vector<double> timeEval{ 1.0, -1.1, +2.2, 3.3 };
	Psym symbol3(name3, timeEval);
	symbol3.print(output);
	BOOST_TEST(!output.is_empty(false));
	stringstream timeEvalString;
	// we compare strings, this is tricky with the length and the single digits
	// workaround. can we do it better. The lexical cast is what the print function does. Is that a good idea?
	// The separator is templated for Psym. SHould we be able to read it?
	timeEvalString << boost::lexical_cast<string>(1.0) << ';' << boost::lexical_cast<string>(-1.1) << ';' << boost::lexical_cast<string>(2.2) << ';' << boost::lexical_cast<string>(3.3);
	string const testString = "name=name3\ntime_eval=" + timeEvalString.str() + "\norder=1\n";
	BOOST_TEST(output.is_equal(testString));


	// construct with name and time evaluation vector and order
	unsigned int const testOrder = 99;
	Psym symbol4(name4, timeEval, testOrder);
	symbol4.print(output);

	BOOST_TEST(!output.is_empty(false));

	string const testString2 = "name=name4\ntime_eval=" + timeEvalString.str() + "\norder=99\n";
	BOOST_TEST(output.is_equal(testString2));
}
