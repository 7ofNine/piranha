#define BOOST_TEST_MODULE vector_key Test
#include "boost/test/included/unit_test.hpp"
#include "boost/test/tools/output_test_stream.hpp"


#include "piranha.h"

using namespace std;
using namespace piranha;
using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace {
	// create a derived type to test the vector_key

	template < class T, int Position>
	class Intermediate : public VectorKey < T, Position, Intermediate<T, Position>>
	{
		using Ancestor = VectorKey < T, Position, Intermediate<T, Position>>;
	public:
		Intermediate() = default;
		Intermediate(Intermediate const &) = default;
		Intermediate(Intermediate &&) = default;

		Intermediate & operator=(Intermediate const &) = default;
		Intermediate & operator=(Intermediate &&) = default;

		template <class ArgsTuple>
		Intermediate(const Psym & p, const int & n, ArgsTuple const & a):Ancestor(p, n, a) {}

		// tests for protected members. SHould we do that?
		std::size_t pelementsHasher() { return this->elementsHasher(); }
		bool pElementsAreZero() const  { return this->elementsAreZero(); }
		void pprintElements(std::ostream & outstream) const { this->printElements(outstream); }
		void pprintElementsSorted(std::ostream &outstream, std::vector<std::pair<bool, std::size_t> > positions) const
		{
			this->printElementsSorted(outstream, positions);
		}

	};

	using KeyType = Intermediate<int, 0>;
	using KeyType1 = Intermediate<int, 1>;
}

void setup() { PsymManager::clear(); }
BOOST_AUTO_TEST_CASE(simple_construction, *utf::fixture(&setup))
{
	KeyType constructed;
	BOOST_TEST((constructed.size() == 0)); // this is also a test for size()

}

BOOST_AUTO_TEST_CASE(copy_construction, *utf::fixture(&setup))
{
	KeyType constructed;
	BOOST_CHECK_NO_THROW(KeyType tmp = KeyType());    
	BOOST_CHECK_NO_THROW(KeyType tmp = KeyType(KeyType()));
	BOOST_CHECK_NO_THROW(KeyType tmp(constructed));
}


BOOST_AUTO_TEST_CASE(construct_with_psym, *utf::fixture(&setup))
{
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

// how to test that tested through derived class TrigVector
// does not exist fro ExpoVector. Is that really needed??
BOOST_AUTO_TEST_CASE(construct_from_different_position, *utf::fixture(&setup))
{
	//KeyType1 temp1;  
	//temp1.resize(3);

	//KeyType  tmp(temp1);
	//BOOST_FAIL("Not implemented");
	
}


BOOST_AUTO_TEST_CASE(resize, *utf::fixture(&setup))
{
	KeyType constructed;
	constructed.resize(5);
	BOOST_TEST((constructed.size() == 5)); // this also a test for size()
}


BOOST_AUTO_TEST_CASE(swapkey, *utf::fixture(&setup))
{
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;

	KeyType temp1;
	
	temp1.swap(temp);
	BOOST_TEST(temp.size() == 0);
	BOOST_TEST(temp1.size() == 3);
	BOOST_TEST(temp1[0] == 1);
	BOOST_TEST(temp1[1] == 2);
	BOOST_TEST(temp1[2] == 3);
}

BOOST_AUTO_TEST_CASE(needs_padding, *utf::fixture(&setup))
{
	Psym t1("t1");
	Psym t2("t2");

	VectorPsym polySym = { t1, t2 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	KeyType temp;
	temp.resize(1);
	BOOST_TEST(temp.needsPadding(polyOnlyArgs));
}

BOOST_AUTO_TEST_CASE(isInsertable, *utf::fixture(&setup))
{
	Psym t1("t1");
	Psym t2("t2");

	VectorPsym polySym = { t1, t2 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	KeyType temp;
	temp.resize(1);
	BOOST_TEST(temp.isInsertable(polyOnlyArgs));
	temp.resize(2);
	BOOST_TEST(temp.isInsertable(polyOnlyArgs));
	temp.resize(3);
	BOOST_TEST(!temp.isInsertable(polyOnlyArgs));
}

BOOST_AUTO_TEST_CASE(atoms, *utf::fixture(&setup))
{
	KeyType temp;
	BOOST_TEST(temp.atoms() ==1 );
}


BOOST_AUTO_TEST_CASE(padRight, *utf::fixture(&setup))
{
	KeyType temp;
	Psym t1("t1");
	Psym t2("t2");

	VectorPsym polySym = { t1, t2 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	BOOST_CHECK_NO_THROW(temp.padRight(polyOnlyArgs));
	BOOST_TEST(temp.size() == 2);

	BOOST_CHECK_NO_THROW(temp.padRight(polyOnlyArgs));// already the same size

	temp.resize(5);
	BOOST_CHECK_THROW(temp.padRight(polyOnlyArgs), assertion_error);
}

BOOST_AUTO_TEST_CASE(applyLayout, *utf::fixture(&setup))
{	// seems only to be used from namedSeries
	// needed but not really used in the method, why is actually in the template
	Psym t1("t1");
	Psym t2("t2");

	VectorPsym polySym = { t1, t2 };

	boost::tuple<VectorPsym> polyOnlyArgs(polySym);

	// the layoutTuple is defined ad template in named_series_def.h
	LayoutTuple<boost::tuple<VectorPsym>>::Type layoutTuple;

	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	KeyType tempSaved(temp);

	// no element at all
	LayoutElement nullPair = std::make_pair(false, 3);
	Layout layout;
	layoutTuple.get<0>() = Layout();
	// layouTuple doesn't have enough elements
	BOOST_CHECK_THROW(temp.applyLayout(layoutTuple, polyOnlyArgs), assertion_error);
	
	temp = tempSaved;
	// no element are selected
	layoutTuple.get<0>() = Layout{ nullPair, nullPair, nullPair };
	temp.applyLayout(layoutTuple, polyOnlyArgs);
	BOOST_TEST(temp[0] == 0);
	BOOST_TEST(temp[1] == 0);
	BOOST_TEST(temp[2] == 0);

	// move element 2 to 0 everything else to 0
	LayoutElement pair2 = std::make_pair(true, 2);
	layoutTuple.get<0>() = Layout{ pair2, nullPair, nullPair };
	temp = tempSaved;
	temp.applyLayout(layoutTuple, polyOnlyArgs);
	BOOST_TEST(temp[0] == 3);
	BOOST_TEST(temp[1] == 0);
	BOOST_TEST(temp[2] == 0);

	// reverse the order of elements
	LayoutElement pair0 = std::make_pair(true, 0);
	LayoutElement pair1 = std::make_pair(true, 1);
	layoutTuple.get<0>() = Layout{ pair2, pair1, pair0 };
	temp = tempSaved;
	temp.applyLayout(layoutTuple, polyOnlyArgs);
	BOOST_TEST(temp[0] == 3);
	BOOST_TEST(temp[1] == 2);
	BOOST_TEST(temp[2] == 1);
}

BOOST_AUTO_TEST_CASE(trimTest, *utf::fixture(&setup))
{
	// trimflags is an Ntuple of vector<bool>,
	// true: don't trim (the key value at that index != 0
	// false : trim, the key value at that index == 0
	// they are basically only used in named_series but get transported through from there to this level
	//typedef typename NTuple<std::vector<bool>, Derived::echelonLevel + 1>::Type TrimFlagsType;
	using TrimFlagsType = NTuple<std::vector<bool>, 1>::Type;
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	std::vector<bool> flags(4, false);
	TrimFlagsType trimFlags;
	trimFlags.get<0>() = flags;

	// different sizes
	BOOST_CHECK_THROW(temp.trimTest(trimFlags), assertion_error);

	std::vector<bool> const flags3(3, false);
	trimFlags.get<0>() = flags3;
	// don't trim any 
	BOOST_CHECK_NO_THROW(temp.trimTest(trimFlags));
	BOOST_TEST(trimFlags.get<0>()[0] == true);
	BOOST_TEST(trimFlags.get<0>()[1] == true);
	BOOST_TEST(trimFlags.get<0>()[2] == true);

	// set one to 0 i.e. make it "trimmable"
	temp[1] = 0;
	trimFlags.get<0>() = flags3;
	// don't trim any 
	BOOST_CHECK_NO_THROW(temp.trimTest(trimFlags));
	BOOST_TEST(trimFlags.get<0>()[0] == true);
	BOOST_TEST(trimFlags.get<0>()[1] == false);
	BOOST_TEST(trimFlags.get<0>()[2] == true);

	// have all trimmed
	temp[0] = 0;
	temp[1] = 0;
	temp[2] = 0;
	trimFlags.get<0>() = flags3;
	// trim them all 
	BOOST_CHECK_NO_THROW(temp.trimTest(trimFlags));
	BOOST_TEST(trimFlags.get<0>()[0] == false);
	BOOST_TEST(trimFlags.get<0>()[1] == false);
	BOOST_TEST(trimFlags.get<0>()[2] == false);
}

BOOST_AUTO_TEST_CASE(trim, *utf::fixture(&setup))
{

	boost::tuple<VectorPsym> polyOnlyArgs; // ArgsTuple is not used for anything, why is it here?

	using TrimFlagsType = NTuple<std::vector<bool>, 1>::Type;
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 0;
	temp[2] = 3;

	TrimFlagsType trimFlags;
	trimFlags.get<0>() = std::vector<bool>(temp.size(), false);
	// create the proper flags, uses another tested function
	temp.trimTest(trimFlags);

	// trim one out
	KeyType trimmed;
	trimmed = temp.trim(trimFlags, polyOnlyArgs);
	BOOST_TEST(trimmed.size() == 2);
	BOOST_TEST(trimmed[0] = 1);
	BOOST_TEST(trimmed[1] = 3);

	// unequal size
	trimFlags.get<0>() = std::vector<bool>(4, false);
	// bigger
	BOOST_CHECK_THROW(KeyType trimmed = temp.trim(trimFlags, polyOnlyArgs), assertion_error);
	//smaller
	trimFlags.get<0>() = std::vector<bool>(2, false);
	BOOST_CHECK_THROW(KeyType trimmed = temp.trim(trimFlags, polyOnlyArgs), assertion_error);
}


BOOST_AUTO_TEST_CASE(invert, *utf::fixture(&setup))
{
	KeyType temp;
	temp.resize(2);
	temp[0] = 1;
	temp[1] = 2;

	BOOST_CHECK_NO_THROW(temp.invertSign());
	BOOST_TEST(temp.size() = 2);
	BOOST_TEST(temp[0] == -1);
	BOOST_TEST(temp[1] == -2);
}

BOOST_AUTO_TEST_CASE(subscript_operator, *utf::fixture(&setup))
{
	// is this a usefull test we already did that many times above
	KeyType temp;
	temp.resize(2);
	temp[0] = 1;
	temp[1] = 2;

	BOOST_CHECK_THROW(temp[5], assertion_error);
	using Stype = KeyType::size_type;
	int const temp0 = temp[0];
	BOOST_TEST(temp0 == 1);
	
	// non const
	int temp1 = temp[1];
	BOOST_TEST(temp1 == 2);
}

BOOST_AUTO_TEST_CASE(iterators, *utf::fixture(&setup))
{
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	KeyType::const_iterator cbiterator = temp.begin();
	KeyType::const_iterator ceiterator = temp.end();
	BOOST_TEST((*cbiterator) == 1);
	BOOST_TEST((*(ceiterator-1)) == 3);

	KeyType::iterator biterator = temp.begin();
	KeyType::iterator eiterator = temp.end();
	BOOST_TEST((*biterator) == 1);
	BOOST_TEST((*(eiterator-1)) == 3);
}

BOOST_AUTO_TEST_CASE(revlex_comparison, *utf::fixture(&setup))
{
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	KeyType temp1(temp);

	// they are equal i.e. fail
	BOOST_TEST(!temp.revlexComparison(temp1)); // (1,2,3) == (1,2,3)
	// all less than temp1 should pass
	temp[2]--;
	BOOST_TEST(temp.revlexComparison(temp1));  // (1,2,2) < (1,2,3)
	// make one larger than temp1
	temp[1]++;
	BOOST_TEST(temp.revlexComparison(temp1)); // (1,3,2) < (1,2,3) 
	temp[1]--;
	BOOST_TEST(temp.revlexComparison(temp1)); // (1,2,2) < (1,2,3) 
	temp[2]++;
	BOOST_TEST(!temp.revlexComparison(temp1));// (1,3,3) > (1,2,3) 
	// do we need more tests?
	// make unequal size
	temp.resize(5);
	BOOST_CHECK_THROW(temp.lexComparison(temp1), assertion_error);
}


BOOST_AUTO_TEST_CASE(lexComparison, *utf::fixture(&setup))
{
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	KeyType temp1(temp);

	// they are equal i.e. fail
	BOOST_TEST(!temp.lexComparison(temp1)); // (1,2,3) == (1,2,3)
	// all less than temp1 should pass
	temp[0]--;
	temp[1]--;
	temp[2]--;
	BOOST_TEST(temp.lexComparison(temp1));  // (0,1,2) < (1,2,3)
	// make one larger than temp1
	temp[1]++; 
	BOOST_TEST(temp.lexComparison(temp1)); //(0,2,2) < (1,2,3) 
	temp[0]++;
	BOOST_TEST(temp.lexComparison(temp1)); //(1,2,2) < (1,2,3) 
	temp[1]++;
	BOOST_TEST(!temp.lexComparison(temp1)); //(1,3,2) > (1,2,3) 
	// do we need more tests?

	// make unequal size
	temp.resize(5);
	BOOST_CHECK_THROW(temp.lexComparison(temp1), assertion_error);
}

BOOST_AUTO_TEST_CASE(operatorEqual, *utf::fixture(&setup))
{
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	KeyType temp1(temp);

	BOOST_TEST(temp == temp1);
	
	temp[0] = 0;
	BOOST_TEST(!(temp == temp1));

	// different sizes (should fail
	temp1.resize(5);
	temp == temp1;
}

BOOST_AUTO_TEST_CASE(elementsEqualTo, *utf::fixture(&setup))
{
	// Hmm. The code is identical to operator==. wht is this good for
	// even used anywhere

	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;
	KeyType temp1(temp);

	BOOST_TEST(temp.elementsEqualTo(temp1));

	temp[0] = 0;
	BOOST_TEST(!(temp.elementsEqualTo(temp1)));

	// different sizes (should fail
	temp1.resize(5);
	temp.elementsEqualTo(temp1);
}

BOOST_AUTO_TEST_CASE(printElements, *utf::fixture(&setup))
{
	// this is a protected member hence we use a tes wrapper
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;

	output_test_stream output;
	decltype(KeyType::separator) const sep = KeyType::separator;
	BOOST_CHECK_NO_THROW(temp.pprintElements(output));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("1;2;3"));

	temp.invertSign();
	BOOST_CHECK_NO_THROW(temp.pprintElements(output));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("-1;-2;-3"));
}


BOOST_AUTO_TEST_CASE(printElementsSorted, *utf::fixture(&setup))
{
	KeyType temp;
	temp.resize(3);
	temp[0] = 1;
	temp[1] = 2;
	temp[2] = 3;

	output_test_stream output;
	std::vector<std::pair<bool, std::size_t>> sort = { {true, 2}, {true, 0}, {true, 1} };
	BOOST_CHECK_NO_THROW(temp.pprintElementsSorted(output, sort));
	BOOST_TEST(!output.is_empty(false));
	BOOST_TEST(output.is_equal("     3      1      2 "));

	//sizes don't agree
	BOOST_CHECK_THROW(temp.pprintElementsSorted(output, std::vector<std::pair<bool, std::size_t>>()), assertion_error);
}

BOOST_AUTO_TEST_CASE(elementsAreZero, *utf::fixture(&setup))
{
	KeyType temp;

	BOOST_TEST(temp.pElementsAreZero());

	temp.resize(3);
	temp[0] = 0;
	temp[1] = 0;
	BOOST_TEST(temp.pElementsAreZero());
	temp[2] = 0;
	BOOST_TEST(temp.pElementsAreZero());
	temp[2]++;
	BOOST_TEST(!temp.pElementsAreZero());
	temp[0]--;
	temp[2]++;
	BOOST_TEST(!temp.pElementsAreZero());
}

BOOST_AUTO_TEST_CASE(elementsHasher, *utf::fixture(&setup))
{
	//// not a real test. just to eecute the algorithm
	KeyType temp;
	BOOST_TEST(temp.pelementsHasher() == 0);

	temp.resize(3);
	temp[0] = -2;
	temp[1] = 3;
	temp[2] = 5;
	// Can we do something better here?
	BOOST_CHECK_NO_THROW(temp.pelementsHasher());
}

