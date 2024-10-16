#define BOOST_TEST_MODULE testexample
#include <boost/test/unit_test.hpp>

int add(int i, int j) { return i+j;}

BOOST_AUTO_TEST_CASE (testexample)
{
  // six ways to detect and report the sme error
  BOOST_CHECK( add(2,2)==4);            // continues on error

  BOOST_REQUIRE( add(2,2)==4);          // throws on error

  if (add(2,2) !=4)
    BOOST_ERROR("Here is an error..."); // continues on error
  
  if (add(2,2) !=4)
    BOOST_FAIL("Here is an error...");  // throws on error

  BOOST_CHECK_MESSAGE( add(2,2) ==4,    // continues on error  
                       "add(...) result: " << add(2,2)); 

  BOOST_CHECK_EQUAL( add(2,2), 4);      // continues on error
    
}


