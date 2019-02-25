#define BOOST_TEST_MODULE BaseSeries Test
#include "boost/test/included/unit_test.hpp"
#include "boost/test/output_test_stream.hpp"


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


// argTuples are introduced using Psym but their connection to names is not actuall being used! That is bussines of named series

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

    using STerm = Monomial<double_cf, ExpoVector<int, 0>, '|', std::allocator<char>>;
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