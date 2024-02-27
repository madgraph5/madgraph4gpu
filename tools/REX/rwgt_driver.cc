//==========================================================================
// Copyright (C) 2023-2024 CERN
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Written by: Z. Wettersten (Jan 2024) for the MG5aMC CUDACPP plugin.
//==========================================================================
//==========================================================================
// This file has been automatically generated for C++ Standalone by
//%(info_lines)s
//==========================================================================
//==========================================================================
// Driver for reweighting events for processes
//%(multiprocess_lines)s
//--------------------------------------------------------------------------

#include "rwgt_instance.h"
#include <cstdlib>
#include <typeinfo>
//%(include_lines)s

int usage( char* argv0, int ret = 1 )
{
    std::cout << "Usage: " << argv0
        << " [--lhefile=\"/YOUR/PATH/HERE\"|-lhe=\"/YOUR/PATH/HERE\"] [--rwgtcard=/YOUR/PATH/HERE|-rwgt=\"/YOUR/PATH/HERE\"]\n"
        << "[--output=/YOUR/PATH/HERE\"|-out=\"/YOUR/PATH/HERE\"]\n" << "[--param_card=/YOUR/PATH/HERE\"|-slha=\"/YOUR/PATH/HERE\"]\n";
    std::cout << "\n";
    std::cout << "The LHE file path should be with respect to the directory you are running\n";
    std::cout << "this program from, and similarly the rwgt_card should be as well.\n";
    return ret;
}


int main( int argc, char** argv ){
    std::string lheFilePath;
    std::string rwgtCardPath;
    std::string outputPath;
    std::string slhaPath;

    if (argc < 2){
        return usage( argv[0] );
    }

    // READ COMMAND LINE ARGUMENTS
    for( int i = 1; i <= argc; i++ )
    {
        auto currArg = std::string( argv[i] );
        if( currArg.substr(0,9) == "--lhefile" || currArg.substr(0,4) == "-lhe" )
        {
            lheFilePath = currArg.substr( currArg.find( "=" ) + 1 ); 
        }
        else if( currArg.substr(0,10) == "--rwgtcard" || currArg.substr(0,5) == "-rwgt" )
        {
            rwgtCardPath = currArg.substr( currArg.find( "=" ) + 1 );
        } else if( currArg.substr(0,8) == "--output" || currArg.substr(0,4) == "-out" ){
            outputPath = currArg.substr( currArg.find( "=" ) + 1 );
        } else if (currArg.substr(0,12) == "--param_card" || currArg.substr(0,5) == "-slha" ){
            slhaPath = currArg.substr( currArg.find( "=" ) + 1 );
        }
        {
            return usage( argv[0] );
        }
    }

    if( lheFilePath.empty() || rwgtCardPath.empty() ){
        return usage( argv[0] );
    }

    std::string currPath = argv[0];

    size_t slashPos = currPath.find_last_of( "/" ); 
    bool onWindows = false;
    if( slashPos == std::string::npos ){ slashPos = currPath.find_last_of( "\\" ); onWindows = true; }
    if( slashPos == std::string::npos )
        throw std::runtime_error( "Failed to determine current working directory -- need to know where program is run from to identify where to pull and push param_card.dat." );

    if( slhaPath.empty() ){
    if( onWindows ){
        if( currPath.substr( currPath.find_last_of("\\", slashPos - 1) + 1, 2 ) == "P1" ){
            slhaPath = "..\\..\\Cards\\param_card.dat";
        } else{
            slhaPath = "\\Cards\\param_card.dat";
        }
    } else {
        if( currPath.substr( currPath.find_last_of("/", slashPos - 1) + 1, 2 ) == "P1" ){
            slhaPath = "../../Cards/param_card.dat";
        } else {
            slhaPath = "/Cards/param_card.dat";
        }
    }}
    

    // ZW : include rwgt_instances(s)
//%(rwgt_runners)s

//    std::vector<rwgt::instance> runSet = {%(run_set)s};
    std::vector<rwgt::instance> runSet;
    REX::teaw::rwgtFiles fileCol( lheFilePath, slhaPath, rwgtCardPath );
    fileCol.initCards();

    REX::teaw::ampCall subProcSet;

    for( auto proc : runSet ){
        subProcSet.insert( REX::teaw::ampPair( proc.procEvent, proc.bridgeCall ) );
    }

    //auto bridgeCont = fbridgeRunner( fileCol.getLhe() );

    //std::function<std::shared_ptr<std::vector<FORTRANFPTYPE>>( std::vector<double>&, std::vector<double>& )> scatteringAmplitude = bridgeCont.scatAmp;
    REX::teaw::rwgtRunner driver( fileCol, subProcSet );


    driver.runRwgt( outputPath ); 
 
    return 0;

}