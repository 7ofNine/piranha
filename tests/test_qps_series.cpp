#define BOOST_TEST_MODULE qps_series_Test
#include "boost/test/included/unit_test.hpp"

// qps_serises_test test the specific qps instantiation of PoissonSeries
// This also allows to test methods in underlying ares that are either in the wrong layer
// originally or basically can not be reached
// als used for experiments of tests

// Should we create a suit eand put all the derived classes in there?
#include "piranha.h"

using namespace piranha;
using namespace piranha::manipulators;

BOOST_AUTO_TEST_CASE(construction_test)
{
    qps testSeries("test_qps_data01.qps");
    auto testresult = testSeries.partial("t1", 1);

}