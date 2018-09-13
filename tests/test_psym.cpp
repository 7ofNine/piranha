#define BOOST_TEST_MODULE Psym Test
#include "boost/test/included/unit_test.hpp"

#include <string>

#include "piranha.h"

using namespace std;
using namespace piranha;

BOOST_AUTO_TEST_CASE(construction_test)
{
	string const name1("name1");
	string const name2("name2");
	int const defaultOrder = 1; // the default order set ib constructor
	
	// construct with name only
	Psym symbol1(name1);
	BOOST_TEST(symbol1.getName() == "name1");
	BOOST_TEST(symbol1.getTimeEval().empty());
	BOOST_TEST(symbol1.order() == defaultOrder);

	// construct with name and order
	Psym symbol2(name2, 3);
	BOOST_TEST(symbol2.getName() == "name2");
	BOOST_TEST(symbol2.getTimeEval().empty());
	BOOST_TEST(symbol2.order() == 3);

	// construct with name and time evaluation vector
	string const name3("name3");
	vector<double> timeEval{ 1.0, -1.1, +2.2, 3.3 };
	// the string deliberately has blanks in the middle, to check string normalization before casting
	string const timeEvalString("1.0; -1.1;+2.2 ; 3.3 "); // can we get to the separator ";" programatically ?
	Psym symbol3(name3, timeEval);
	BOOST_TEST(symbol3.getName() == name3);
	BOOST_TEST(symbol3.getTimeEval() == timeEval);
	BOOST_TEST(symbol3.order() = defaultOrder);

	// construct with name and time evaluation vector and order
	string const name4("name4");
	unsigned int testOrder = 99;
	Psym symbol4(name4, timeEval, testOrder);
	BOOST_TEST(symbol4.getName() == name4);
	BOOST_TEST(symbol4.getTimeEval() == timeEval);
	BOOST_TEST(symbol4.order() == testOrder);

	//construct with name and time evaluation string
	string const name5("name5");
	Psym symbol5(name5, timeEvalString);
	BOOST_TEST(symbol5.getName() == name5);
	BOOST_TEST(symbol5.getTimeEval() == timeEval);
	BOOST_TEST(symbol5.order() == defaultOrder);


	//construct with name and time evaluation string and order
	string const name6("name6");
	Psym symbol6(name6, timeEvalString, testOrder);
	BOOST_TEST(symbol6.getName() == name6);
	BOOST_TEST(symbol6.getTimeEval() == timeEval);
	BOOST_TEST(symbol6.order() == testOrder);

	//construct with name and empty time evaluation string and order
	string const name7("name7");
	string const emptyTimeEvalString("");
	Psym symbol7(name7, emptyTimeEvalString, testOrder);
	BOOST_TEST(symbol7.getName() == name7);
	BOOST_TEST(symbol7.getTimeEval().empty());
	BOOST_TEST(symbol7.order() == testOrder);

}