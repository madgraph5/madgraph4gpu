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
#ifndef _LIBCOMP_
#define _LIBCOMP_
#endif
#include "rwgt_instance.h"
#include "fbridge.cc"

// ZW: SET UP NAMESPACE
namespace %(process_namespace)s{
//namespace dummy{

    std::shared_ptr<std::vector<FORTRANFPTYPE>> amp( int& nEvt, int& nPar, int& nMom, std::vector<FORTRANFPTYPE>& momenta, std::vector<FORTRANFPTYPE>& alphaS, std::vector<FORTRANFPTYPE>& rndHel, std::vector<FORTRANFPTYPE>& rndCol, std::vector<int>& selHel, std::vector<int>& selCol, unsigned int& chanId ){
        CppObjectInFortran *bridgeInst;
        auto evalScatAmps = std::make_shared<std::vector<FORTRANFPTYPE>>( nEvt );
        fbridgecreate_( &bridgeInst, &nEvt, &nPar, &nMom );
        fbridgesequence_( &bridgeInst, &momenta.at(0), &alphaS.at(0), &rndHel[0], &rndCol[0], &chanId, &evalScatAmps->at(0), &selHel[0], &selCol[0] );
        fbridgedelete_( &bridgeInst );
        return evalScatAmps;
    }

    rwgt::fBridge bridgeConstr( std::vector<REX::event>& process, unsigned int warpSize = 32 ){
        rwgt::fBridge constrBridge =  rwgt::fBridge( process, warpSize );
        rwgt::bridgeWrapper amplitude = amp; 
        constrBridge.setBridge( amplitude );
        return constrBridge;
    }

    rwgt::fBridge bridgeConstr(){
        rwgt::fBridge constrBridge =  rwgt::fBridge();
        rwgt::bridgeWrapper amplitude = amp; 
        constrBridge.setBridge( amplitude );
        return constrBridge;
    }

    std::shared_ptr<std::vector<size_t>> procSort( std::string_view status, std::vector<std::string_view> arguments ){
        std::vector<std::vector<std::string_view>> initPrts = {%(init_prt_ids)s};
        std::vector<std::vector<std::string_view>> finPrts = {%(fin_prt_ids)s};
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

    REX::eventSet eventSetConstr( std::vector<REX::event>& process ){
        REX::eventSet constrSet = REX::eventSet( process );
        REX::eventSetComp compar = checkProc;
        constrSet.setComp( compar );
        return constrSet;
    }

    REX::eventSet getEventSet(){
        std::vector<std::vector<std::pair<std::string,std::string>>> eventVec = {%(process_events)s};
        std::vector<REX::event> process;
        for( auto ev : eventVec ){
            process.push_back( REX::event( ev ) );
        }
        return eventSetConstr( process );
    }

}