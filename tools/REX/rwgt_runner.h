//==========================================================================
// Copyright (C) 2023-2024 CERN
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Written by: Z. Wettersten (June 2024) for the MG5aMC CUDACPP plugin.
//==========================================================================
//==========================================================================
// This file has been automatically generated for the CUDACPP plugin by
%(info_lines)s
//==========================================================================
//==========================================================================
// A class for reweighting matrix elements for
%(process_lines)s
//--------------------------------------------------------------------------

#ifndef _%(process_namespace)s_RUNNER_H_
#define _%(process_namespace)s_RUNNER_H_

#include "rwgt_instance.h"

namespace %(process_namespace)s {

    std::shared_ptr<std::vector<FORTRANFPTYPE>> amp( int& nEvt, int& nPar, int& nMom, std::vector<FORTRANFPTYPE>& momenta, std::vector<FORTRANFPTYPE>& alphaS, std::vector<FORTRANFPTYPE>& rndHel, std::vector<FORTRANFPTYPE>& rndCol, std::vector<int>& selHel, std::vector<int>& selCol, int& chanId );
    rwgt::fBridge bridgeConstr( std::vector<REX::event>& process, unsigned int warpSize );
    rwgt::fBridge bridgeConstr();
    std::shared_ptr<std::vector<size_t>> procSort( std::string_view status, std::vector<std::string_view> arguments, size_t index = REX::npos );
    bool checkProc( REX::event& process, std::vector<std::string>& relStats );
    REX::eventSet eventSetConstruct( std::vector<REX::event>& process );
    REX::eventSet getEventSet();

}



#endif