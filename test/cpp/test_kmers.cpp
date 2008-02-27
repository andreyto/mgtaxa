#include "mgtaxa/kmers.hpp"
#include "mgtaxa/version.hpp"
#include <boost/test/minimal.hpp>

namespace {

bool test_KmerCounter() {
	MGT::KmerCounter counter(2);
	return true;
}

}

int add( int i, int j ) { return i+j; }

int test_main( int, char *[] )             // note the name!
{
	BOOST_CHECK(test_KmerCounter());
    // six ways to detect and report the same error:
/*    BOOST_CHECK( add( 2,2 ) == 4 );        // #1 continues on error
    BOOST_REQUIRE( add( 2,2 ) == 5 );      // #2 throws on error
    if( add( 2,2 ) != 4 )
      BOOST_ERROR( "Ouch..." );            // #3 continues on error
    if( add( 2,2 ) != 4 )
      BOOST_FAIL( "Ouch..." );             // #4 throws on error
    if( add( 2,2 ) != 4 ) throw "Oops..."; // #5 throws on error

    return add( 2, 2 ) == 4 ? 0 : 1;       // #6 returns error code*/
}
