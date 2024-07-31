//==========================================================================
// Copyright (C) 2023-2024 CERN
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Written by: Z. Wettersten (Jan 2024) for the MG5aMC CUDACPP plugin.
//==========================================================================
//==========================================================================
// This file has been automatically generated for the CUDACPP plugin by
%(info_lines)s
//==========================================================================
//==========================================================================
// A class for reweighting matrix elements for
%(process_lines)s
//--------------------------------------------------------------------------

#include "rwgt_instance.h"
#include "fbridge.cc"

// ZW: SET UP NAMESPACE
namespace %(process_namespace)s{
//namespace dummy{

    std::shared_ptr<std::vector<FORTRANFPTYPE>> amp( int& nEvt, int& nPar, int& nMom, std::vector<FORTRANFPTYPE>& momenta, std::vector<FORTRANFPTYPE>& alphaS, std::vector<FORTRANFPTYPE>& rndHel, std::vector<FORTRANFPTYPE>& rndCol, std::vector<int>& selHel, std::vector<int>& selCol, int& chanId ){
        CppObjectInFortran *bridgeInst;
        auto evalScatAmps = std::make_shared<std::vector<FORTRANFPTYPE>>( nEvt );
        fbridgecreate_( &bridgeInst, &nEvt, &nPar, &nMom );
        fbridgesequence_( &bridgeInst, &momenta.at(0), &alphaS.at(0), &rndHel[0], &rndCol[0], &chanId, &evalScatAmps->at(0), &selHel[0], &selCol[0] );
        fbridgedelete_( &bridgeInst );
        return evalScatAmps;
    }

    rwgt::fBridge& bridgeConstr( std::vector<REX::event>& process, unsigned int warpSize = 32 ){
        rwgt::fBridge constrBridge =  rwgt::fBridge( process, warpSize );
        rwgt::bridgeWrapper amplitude = amp; 
        constrBridge.setBridge( amplitude );
        return constrBridge;
    }

    std::shared_ptr<std::vector<size_t>> procSort( std::string_view status, std::vector<std::string_view> arguments ){
        std::vector<std::vector<std::string>> initPrts = {%(init_prt_ids)s};
        std::vector<std::vector<std::string>> finPrts = {%(fin_prt_ids)s};
//        std::vector<std::string_view> initPrts = {"-1"};
//        std::vector<std::string_view> finPrts = {"1"};
        std::shared_ptr<std::vector<size_t>> refOrder;
        if( status == "-1" ){
            for( auto& prts : initPrts ){
                refOrder = REX::getRefOrder( prts, arguments );
                if( refOrder->at(refOrder->size() - 1) != REX::npos ){ break; }
            }
            return refOrder;
        }
        else if( status == "1" ){
            for( auto& prts : finPrts ){
                refOrder = REX::getRefOrder( prts, arguments );
                if( refOrder->at(refOrder->size() - 1) != REX::npos ){ break; }
            }
            return refOrder;
        }
        return REX::stoiSort( arguments );
    }

    bool checkProc( REX::event& process, std::vector<std::string>& relStats ){
        REX::statSort locSort = procSort;
        auto order = process.getProcOrder( locSort );
        for( auto stat : relStats ){
            auto currPts = order.at( stat );
            if( currPts[currPts.size() - 1 ] == REX::npos ){ return false; }
        }
        return true;
    }

    REX::eventSet& eventSetConstr( std::vector<REX::event>& process ){
        REX::eventSet constrSet = REX::eventSet( process );
        REX::eventSetComp compar = checkProc;
        constrSet.setComp( compar );
        return constrSet;
    }

    REX::eventSet& getEventSet(){
        std::vector<std::vector<std::pair<std::string,std::string>>> eventVec = {%(process_events)s};
        std::vector<REX::event> process;
        for( auto& ev : eventVec ){
            process.push_back( REX::event( ev ) );
        }
        return eventSetConstr( process );
    }

    struct fbridgeRunner{
        std::vector<FORTRANFPTYPE> rndHel;
        std::vector<FORTRANFPTYPE> rndCol;
        std::vector<int> selHel;
        std::vector<int> selCol;
        CppObjectInFortran *fBridge;
        const unsigned int chanId = 0;
        const int nMom = 4;
        int nWarpRemain;
        int nEvt;
        int fauxNEvt;
        int nPar;
        bool setup = false;
        fbridgeRunner(){}
        fbridgeRunner( REX::event& process ){
            nPar = process.getPrts().size();
        }
        void runnerSetup( unsigned int& noEvts, unsigned int warpSize = 32){
            if( setup ){ return; }
            nEvt = noEvts;
            nWarpRemain = rwgt::warpRemain( nEvt, warpSize );
            fauxNEvt = nEvt + nWarpRemain;
            rndHel = std::vector<FORTRANFPTYPE>( fauxNEvt, 0. );
            rndCol = std::vector<FORTRANFPTYPE>( fauxNEvt, 0. );
            selHel = std::vector<int>( fauxNEvt, 0 );
            selCol = std::vector<int>( fauxNEvt, 0 );
            setup = true;
        }
        void runnerSetup( std::vector<FORTRANFPTYPE>& evVec, unsigned int warpSize = 32){
            if( setup ){ return; }
            nEvt = evVec.size();
            nWarpRemain = rwgt::warpRemain( nEvt, warpSize );
            fauxNEvt = nEvt + nWarpRemain;
            rndHel = std::vector<FORTRANFPTYPE>( fauxNEvt, 0. );
            rndCol = std::vector<FORTRANFPTYPE>( fauxNEvt, 0. );
            selHel = std::vector<int>( fauxNEvt, 0 );
            selCol = std::vector<int>( fauxNEvt, 0 );
            setup = true;
        }
        void runnerSetup( std::shared_ptr<std::vector<FORTRANFPTYPE>> evVec, unsigned int warpSize = 32){
            if( setup ){ return; }
            runnerSetup( *evVec, warpSize );
        }
        std::shared_ptr<std::vector<FORTRANFPTYPE>> scatAmp( std::vector<FORTRANFPTYPE>& momenta, std::vector<FORTRANFPTYPE>& alphaS ){
            runnerSetup( alphaS );
            for( auto j = 0 ; j < nWarpRemain ; ++j ){
                alphaS.push_back( 0. );
                for( auto k = 0 ; k < nMom * nPar ; ++k ){
                    momenta.push_back( 0. );
                }
            }
            auto evalScatAmps = std::make_shared<std::vector<FORTRANFPTYPE>>( fauxNEvt );
            fbridgecreate_( &fBridge, &fauxNEvt, &nPar, &nMom );
            fbridgesequence_( &fBridge, &momenta.at(0), &alphaS.at(0), &rndHel[0], &rndCol[0], &chanId, &evalScatAmps->at(0), &selHel[0], &selCol[0] );
            fbridgedelete_( &fBridge );
            alphaS.resize( nEvt );
            momenta.resize( nEvt * nPar * nMom );
            evalScatAmps->resize( nEvt );
            return evalScatAmps;
        }
        std::shared_ptr<std::vector<FORTRANFPTYPE>> scatAmp( std::shared_ptr<std::vector<FORTRANFPTYPE>> momenta, std::shared_ptr<std::vector<FORTRANFPTYPE>> alphaS ){
            return scatAmp( *momenta, *alphaS );
        }
#if defined MGONGPU_FPTYPE_FLOAT
        std::shared_ptr<std::vector<FORTRANFPTYPE>> scatAmp( std::vector<double>& momenta, std::vector<double>& alphaS ){
            auto nuMom = std::vector<float>( nEvt );
            auto nuAlphaS = std::vector<float>( nEvt );
            std::transform( momenta.begin(), momenta.end(), nuMom.begin(), [](double mom){ return static_cast<float>(mom); })
            std::transform( alphaS.begin(), alphaS.end(), nuAlphaS.begin(), [](double gs){ return static_cast<float>(gs); });
            return scatAmp( nuMom, nuAlphaS );
        }
#endif
    }; 

    std::shared_ptr<std::vector<size_t>> thisProcSort( std::string_view& status, std::vector<std::string_view>& arguments ){
        std::vector<std::string_view> initPrts = %(init_prt_ids)s;
        std::vector<std::string_view> finPrts = %(fin_prt_ids)s;
//        std::vector<std::string_view> initPrts = {"-1"};
//        std::vector<std::string_view> finPrts = {"1"};
        if( status == "-1" ){
            return REX::getRefOrder( initPrts, arguments );
        }
        else if( status == "1" ){
            return REX::getRefOrder( finPrts, arguments );
        }
        return REX::stoiSort( arguments );
    }

// ZW: SET UP INPUT LHE BLOCK
// ZW: SET UP REX::event FROM LHE BLOCK
//    auto procEvent = REX::event( procEvent );
//    REX::statSort currProcSort = []( std::string_view stat, std::vector<std::string_view> vec ){ return thisProcSort( stat, vec ); };

    std::vector<std::pair<std::string,std::string>> eventVec = {%(process_event)s};
    REX::event locEv = REX::event( eventVec );
    fbridgeRunner fBridge = fbridgeRunner( locEv );

    REX::teaw::amplitude scatteringAmp = []( std::vector<FORTRANFPTYPE>& momenta, std::vector<FORTRANFPTYPE>& alphaS ){
        return fBridge.scatAmp( momenta, alphaS );
    };

    REX::statSort currProcSort = []( std::string_view stat, std::vector<std::string_view> vec ){ return thisProcSort( stat, vec ); };

    auto runner = rwgt::instance(eventVec, scatteringAmp);
    auto thisProc = runner.process.getProc( currProcSort );
// ZW: SET UP WRAPPER FOR FORTRAN_BRIDGE

// ZW: SET UP EVALUATION OF MATRIX ELEMENTS FUNCTION


}