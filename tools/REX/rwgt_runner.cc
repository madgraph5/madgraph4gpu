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
#include "fbridge.h"

// ZW: SET UP NAMESPACE
namespace %(process_namespace)s{
//namespace dummy{

    std::vector<std::vector<std::string_view>> getInitPrts(){
        static std::vector<std::vector<std::string_view>> initPrts = {%(init_prt_ids)s};
        return initPrts;
    }

    std::vector<std::vector<std::string_view>> getFinPrts(){
        static std::vector<std::vector<std::string_view>> finPrts = {%(fin_prt_ids)s};
        return finPrts;
    }

    std::shared_ptr<std::vector<FORTRANFPTYPE>> amp( int& nEvt, int& nPar, int& nMom, std::vector<FORTRANFPTYPE>& momenta, std::vector<FORTRANFPTYPE>& alphaS, std::vector<FORTRANFPTYPE>& rndHel, std::vector<FORTRANFPTYPE>& rndCol, std::vector<int>& selHel, std::vector<int>& selCol, unsigned int& chanId, bool& goodHel ){
        CppObjectInFortran *bridgeInst;
        auto evalScatAmps = std::make_shared<std::vector<FORTRANFPTYPE>>( nEvt );
        fbridgecreate_( &bridgeInst, &nEvt, &nPar, &nMom );
        fbridgesequence_nomultichannel_( &bridgeInst, &momenta.at(0), &alphaS.at(0), &rndHel[0], &rndCol[0], &evalScatAmps->at(0), &selHel[0], &selCol[0], &goodHel );
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

    std::shared_ptr<std::vector<size_t>> procSort( std::string_view status, std::vector<std::string_view> arguments, size_t index ){
        std::vector<std::vector<std::string_view>> initPrts = getInitPrts();
        std::vector<std::vector<std::string_view>> finPrts = getFinPrts();
        std::shared_ptr<std::vector<size_t>> refOrder;
    if( index == REX::npos ){
	if( status == "-1" ){
            for( auto& prts : initPrts ){
                refOrder = REX::getRefOrder( prts, arguments );
                if( std::find(refOrder->begin(), refOrder->end(), REX::npos) == refOrder->end() ){ break; }
            }
            return refOrder;
        }
        else if( status == "1" ){
            for( auto& prts : finPrts ){
                refOrder = REX::getRefOrder( prts, arguments );
                if( std::find(refOrder->begin(), refOrder->end(), REX::npos) == refOrder->end() ){ break; }
            }
            return refOrder;
        }
        return REX::stoiSort( arguments );
    }
    else{
        if( index >= initPrts.size() || index >= finPrts.size() ) throw std::runtime_error( "procSort called for out-of-bounds event." );
        if( status == "-1" ){
            refOrder = REX::getRefOrder( initPrts.at(index), arguments );
            return refOrder;
        }
        else if( status == "1" ){
            refOrder = REX::getRefOrder( finPrts.at(index), arguments );
            return refOrder;
        }
        return REX::stoiSort( arguments );
    }
    }

    bool checkProc( REX::event& process, std::vector<std::string>& relStats ){
        size_t no_evts = %(no_events)s;
        auto finPrts = getFinPrts();
        for( size_t k = 0 ; k < no_evts ; ++k ){
            REX::statSort locSort = [ind = k](std::string_view status, std::vector<std::string_view> arguments){
                return procSort( status, arguments, ind );
            };
            auto order = process.getProcOrder( locSort );
            if( order.at("1").size() != finPrts[k].size() ){ continue; }
            for( size_t j = 0 ; j < relStats.size() ; ++j ){
                auto currPts = order.at( relStats[j] );
                if( std::find(currPts.begin(), currPts.end(), REX::npos) != currPts.end() ){ break; }
                if( j == relStats.size() - 1 ){ return true; }
            }
        }
        return false;
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
