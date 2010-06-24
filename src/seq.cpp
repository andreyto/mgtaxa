//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#include "AS_SCM_seqUtils.hpp"
#include "AS_SCM_utils.hpp"
#include "AS_SCM_debug.hpp"

#include <boost/lexical_cast.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

const char progName[] = "phySeq";


namespace PHY
{
	

	
int actionLoad(const po::variables_map& vm)
{
	
	std::string gkpStorePath;
	
    if (vm.count("gkp-store"))
    {
        gkpStorePath = vm["gkp-store"].as<std::string>();
    }
    else
    {
        std::cerr << "Gatekeeper store path was not set.\n";
        return 1;
    }
    clearRangeUpcast(gkpStorePath);
	return 0;
	    
}
} // namespace PHY

int
main(int argc, char **argv)
{

    using namespace AS_SCM;

    // Declare the supported options.
    po::options_description desc("Service utilities for SCM pipeline");
    desc.add_options()
    ("help", "produce help message")
    ("action,a",
     po::value<std::string>(),
     "What action to perform:\n" 
     "clr-cast - set all writable AS_READ_CLEAR_XXX values to the current AS_READ_CLEAR_LATEST")
    ("gkp-store,g",
     po::value<std::string>(),
     "Gatekeeper store path.")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cerr << desc << "\n";
        return 1;
    }

    std::string action;

    if (vm.count("action"))
    {
        action = vm["action"].as<std::string>();
    }
    else
    {
        std::cerr << "Action was not set.\n";
        return 1;
    }
    
    if( action == "clr-cast" ) 
    {
    	return actionClearRangeCast(vm);
    }
    else
    {
    	SCM_LOG << "Unknown 'action' parameter: " << action << "\n";
    	return 1;
    }
	
	return 0;
}
