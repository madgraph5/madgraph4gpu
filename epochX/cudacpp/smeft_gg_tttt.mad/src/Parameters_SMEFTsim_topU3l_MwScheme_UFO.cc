// Copyright (C) 2010 The MadGraph5_aMC@NLO development team and contributors.
// Created by: J. Alwall (Oct 2010) for the MG5aMC CPP backend.
//==========================================================================
// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Modified by: A. Valassi (Sep 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: J. Teig, A. Valassi (2021-2024) for the MG5aMC CUDACPP plugin.
//==========================================================================
// This file has been automatically generated for CUDA/C++ standalone by
// MadGraph5_aMC@NLO v. 3.6.0, 2024-09-30
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "Parameters_SMEFTsim_topU3l_MwScheme_UFO.h"

#include <iomanip>
#include <iostream>

#ifdef MGONGPUCPP_GPUIMPL
using namespace mg5amcGpu;
#else
using namespace mg5amcCpu;
#endif

#ifndef MGONGPU_HARDCODE_PARAM

// Initialize static instance
Parameters_SMEFTsim_topU3l_MwScheme_UFO* Parameters_SMEFTsim_topU3l_MwScheme_UFO::instance = 0;

// Function to get static instance - only one instance per program
Parameters_SMEFTsim_topU3l_MwScheme_UFO*
Parameters_SMEFTsim_topU3l_MwScheme_UFO::getInstance()
{
  if( instance == 0 )
    instance = new Parameters_SMEFTsim_topU3l_MwScheme_UFO();
  return instance;
}

void
Parameters_SMEFTsim_topU3l_MwScheme_UFO::setIndependentParameters( SLHAReader& slha )
{
  zero = 0;                         // define "zero"
  ZERO = 0;                         // define "zero"
  std::vector<int> indices( 2, 0 ); // prepare a vector for indices
  mdl_WH = slha.get_block_entry( "decay", 25, 4.070000e-03 );
  mdl_WW = slha.get_block_entry( "decay", 24, 2.085000e+00 );
  mdl_WZ = slha.get_block_entry( "decay", 23, 2.495200e+00 );
  mdl_WT = slha.get_block_entry( "decay", 6, 1.330000e+00 );
  mdl_ymtau = slha.get_block_entry( "yukawa", 15, 1.777000e+00 );
  mdl_ymm = slha.get_block_entry( "yukawa", 13, 1.056600e-01 );
  mdl_yme = slha.get_block_entry( "yukawa", 11, 5.110000e-04 );
  mdl_ymt = slha.get_block_entry( "yukawa", 6, 1.727600e+02 );
  mdl_ymb = slha.get_block_entry( "yukawa", 5, 4.180000e+00 );
  mdl_ymc = slha.get_block_entry( "yukawa", 4, 1.270000e+00 );
  mdl_yms = slha.get_block_entry( "yukawa", 3, 9.300000e-02 );
  mdl_ymup = slha.get_block_entry( "yukawa", 2, 2.160000e-03 );
  mdl_ymdo = slha.get_block_entry( "yukawa", 1, 4.670000e-03 );
  mdl_linearPropCorrections = slha.get_block_entry( "switches", 1, 0.000000e+00 );
  //aS = slha.get_block_entry( "sminputs", 3, 1.179000e-01 ); // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  mdl_Gf = slha.get_block_entry( "sminputs", 2, 1.166379e-05 );
  mdl_MW = slha.get_block_entry( "sminputs", 1, 8.038700e+01 );
  mdl_LambdaSMEFT = slha.get_block_entry( "smeftcutoff", 1, 1.000000e+03 );
  mdl_cleQt3Im = slha.get_block_entry( "smeftcpv", 53, 0.000000e+00 );
  mdl_cleQt1Im = slha.get_block_entry( "smeftcpv", 52, 0.000000e+00 );
  mdl_cleju3Im = slha.get_block_entry( "smeftcpv", 51, 0.000000e+00 );
  mdl_cleju1Im = slha.get_block_entry( "smeftcpv", 50, 0.000000e+00 );
  mdl_clebQIm = slha.get_block_entry( "smeftcpv", 49, 0.000000e+00 );
  mdl_cledjIm = slha.get_block_entry( "smeftcpv", 48, 0.000000e+00 );
  mdl_ceBIm = slha.get_block_entry( "smeftcpv", 47, 0.000000e+00 );
  mdl_ceWIm = slha.get_block_entry( "smeftcpv", 46, 0.000000e+00 );
  mdl_ceHIm = slha.get_block_entry( "smeftcpv", 45, 0.000000e+00 );
  mdl_cQtQb8Im = slha.get_block_entry( "smeftcpv", 44, 0.000000e+00 );
  mdl_cQtQb1Im = slha.get_block_entry( "smeftcpv", 43, 0.000000e+00 );
  mdl_cjtQd8Im = slha.get_block_entry( "smeftcpv", 42, 0.000000e+00 );
  mdl_cjtQd1Im = slha.get_block_entry( "smeftcpv", 41, 0.000000e+00 );
  mdl_cQujb8Im = slha.get_block_entry( "smeftcpv", 40, 0.000000e+00 );
  mdl_cQujb1Im = slha.get_block_entry( "smeftcpv", 39, 0.000000e+00 );
  mdl_cjuQb8Im = slha.get_block_entry( "smeftcpv", 38, 0.000000e+00 );
  mdl_cjuQb1Im = slha.get_block_entry( "smeftcpv", 37, 0.000000e+00 );
  mdl_cQtjd8Im = slha.get_block_entry( "smeftcpv", 36, 0.000000e+00 );
  mdl_cQtjd1Im = slha.get_block_entry( "smeftcpv", 35, 0.000000e+00 );
  mdl_cjujd81Im = slha.get_block_entry( "smeftcpv", 34, 0.000000e+00 );
  mdl_cjujd11Im = slha.get_block_entry( "smeftcpv", 33, 0.000000e+00 );
  mdl_cjujd8Im = slha.get_block_entry( "smeftcpv", 32, 0.000000e+00 );
  mdl_cjujd1Im = slha.get_block_entry( "smeftcpv", 31, 0.000000e+00 );
  mdl_cjQbd8Im = slha.get_block_entry( "smeftcpv", 30, 0.000000e+00 );
  mdl_cjQbd1Im = slha.get_block_entry( "smeftcpv", 29, 0.000000e+00 );
  mdl_cjQtu8Im = slha.get_block_entry( "smeftcpv", 28, 0.000000e+00 );
  mdl_cjQtu1Im = slha.get_block_entry( "smeftcpv", 27, 0.000000e+00 );
  mdl_cutbd8Im = slha.get_block_entry( "smeftcpv", 26, 0.000000e+00 );
  mdl_cutbd1Im = slha.get_block_entry( "smeftcpv", 25, 0.000000e+00 );
  mdl_cHtbIm = slha.get_block_entry( "smeftcpv", 24, 0.000000e+00 );
  mdl_cHudIm = slha.get_block_entry( "smeftcpv", 23, 0.000000e+00 );
  mdl_cbHIm = slha.get_block_entry( "smeftcpv", 22, 0.000000e+00 );
  mdl_cdHIm = slha.get_block_entry( "smeftcpv", 21, 0.000000e+00 );
  mdl_ctHIm = slha.get_block_entry( "smeftcpv", 20, 0.000000e+00 );
  mdl_cuHIm = slha.get_block_entry( "smeftcpv", 19, 0.000000e+00 );
  mdl_cbBIm = slha.get_block_entry( "smeftcpv", 18, 0.000000e+00 );
  mdl_cdBIm = slha.get_block_entry( "smeftcpv", 17, 0.000000e+00 );
  mdl_cbWIm = slha.get_block_entry( "smeftcpv", 16, 0.000000e+00 );
  mdl_cdWIm = slha.get_block_entry( "smeftcpv", 15, 0.000000e+00 );
  mdl_cbGIm = slha.get_block_entry( "smeftcpv", 14, 0.000000e+00 );
  mdl_cdGIm = slha.get_block_entry( "smeftcpv", 13, 0.000000e+00 );
  mdl_ctBIm = slha.get_block_entry( "smeftcpv", 12, 0.000000e+00 );
  mdl_cuBIm = slha.get_block_entry( "smeftcpv", 11, 0.000000e+00 );
  mdl_ctWIm = slha.get_block_entry( "smeftcpv", 10, 0.000000e+00 );
  mdl_cuWIm = slha.get_block_entry( "smeftcpv", 9, 0.000000e+00 );
  mdl_ctGIm = slha.get_block_entry( "smeftcpv", 8, 0.000000e+00 );
  mdl_cuGIm = slha.get_block_entry( "smeftcpv", 7, 0.000000e+00 );
  mdl_cHWBtil = slha.get_block_entry( "smeftcpv", 6, 0.000000e+00 );
  mdl_cHBtil = slha.get_block_entry( "smeftcpv", 5, 0.000000e+00 );
  mdl_cHWtil = slha.get_block_entry( "smeftcpv", 4, 0.000000e+00 );
  mdl_cHGtil = slha.get_block_entry( "smeftcpv", 3, 0.000000e+00 );
  mdl_cWtil = slha.get_block_entry( "smeftcpv", 2, 0.000000e+00 );
  mdl_cGtil = slha.get_block_entry( "smeftcpv", 1, 0.000000e+00 );
  mdl_cleQt3Re = slha.get_block_entry( "smeft", 129, 0.000000e+00 );
  mdl_cleju3Re = slha.get_block_entry( "smeft", 128, 0.000000e+00 );
  mdl_cleQt1Re = slha.get_block_entry( "smeft", 127, 0.000000e+00 );
  mdl_cleju1Re = slha.get_block_entry( "smeft", 126, 0.000000e+00 );
  mdl_clebQRe = slha.get_block_entry( "smeft", 125, 0.000000e+00 );
  mdl_cledjRe = slha.get_block_entry( "smeft", 124, 0.000000e+00 );
  mdl_cle = slha.get_block_entry( "smeft", 123, 0.000000e+00 );
  mdl_cbl = slha.get_block_entry( "smeft", 122, 0.000000e+00 );
  mdl_cld = slha.get_block_entry( "smeft", 121, 0.000000e+00 );
  mdl_ctl = slha.get_block_entry( "smeft", 120, 0.000000e+00 );
  mdl_clu = slha.get_block_entry( "smeft", 119, 0.000000e+00 );
  mdl_cQe = slha.get_block_entry( "smeft", 118, 0.000000e+00 );
  mdl_cje = slha.get_block_entry( "smeft", 117, 0.000000e+00 );
  mdl_cbe = slha.get_block_entry( "smeft", 116, 0.000000e+00 );
  mdl_ced = slha.get_block_entry( "smeft", 115, 0.000000e+00 );
  mdl_cte = slha.get_block_entry( "smeft", 114, 0.000000e+00 );
  mdl_ceu = slha.get_block_entry( "smeft", 113, 0.000000e+00 );
  mdl_cee = slha.get_block_entry( "smeft", 112, 0.000000e+00 );
  mdl_cQl3 = slha.get_block_entry( "smeft", 111, 0.000000e+00 );
  mdl_cQl1 = slha.get_block_entry( "smeft", 110, 0.000000e+00 );
  mdl_clj3 = slha.get_block_entry( "smeft", 109, 0.000000e+00 );
  mdl_clj1 = slha.get_block_entry( "smeft", 108, 0.000000e+00 );
  mdl_cll1 = slha.get_block_entry( "smeft", 107, 0.000000e+00 );
  mdl_cll = slha.get_block_entry( "smeft", 106, 0.000000e+00 );
  mdl_cHe = slha.get_block_entry( "smeft", 105, 0.000000e+00 );
  mdl_cHl3 = slha.get_block_entry( "smeft", 104, 0.000000e+00 );
  mdl_cHl1 = slha.get_block_entry( "smeft", 103, 0.000000e+00 );
  mdl_ceBRe = slha.get_block_entry( "smeft", 102, 0.000000e+00 );
  mdl_ceWRe = slha.get_block_entry( "smeft", 101, 0.000000e+00 );
  mdl_ceHRe = slha.get_block_entry( "smeft", 100, 0.000000e+00 );
  mdl_cQtQb8Re = slha.get_block_entry( "smeft", 99, 0.000000e+00 );
  mdl_cQtQb1Re = slha.get_block_entry( "smeft", 98, 0.000000e+00 );
  mdl_cjtQd8Re = slha.get_block_entry( "smeft", 97, 0.000000e+00 );
  mdl_cjtQd1Re = slha.get_block_entry( "smeft", 96, 0.000000e+00 );
  mdl_cQujb8Re = slha.get_block_entry( "smeft", 95, 0.000000e+00 );
  mdl_cQujb1Re = slha.get_block_entry( "smeft", 94, 0.000000e+00 );
  mdl_cjuQb8Re = slha.get_block_entry( "smeft", 93, 0.000000e+00 );
  mdl_cjuQb1Re = slha.get_block_entry( "smeft", 92, 0.000000e+00 );
  mdl_cQtjd8Re = slha.get_block_entry( "smeft", 91, 0.000000e+00 );
  mdl_cQtjd1Re = slha.get_block_entry( "smeft", 90, 0.000000e+00 );
  mdl_cjujd81Re = slha.get_block_entry( "smeft", 89, 0.000000e+00 );
  mdl_cjujd11Re = slha.get_block_entry( "smeft", 88, 0.000000e+00 );
  mdl_cjujd8Re = slha.get_block_entry( "smeft", 87, 0.000000e+00 );
  mdl_cjujd1Re = slha.get_block_entry( "smeft", 86, 0.000000e+00 );
  mdl_cjQbd8Re = slha.get_block_entry( "smeft", 85, 0.000000e+00 );
  mdl_cjQbd1Re = slha.get_block_entry( "smeft", 84, 0.000000e+00 );
  mdl_cjQtu8Re = slha.get_block_entry( "smeft", 83, 0.000000e+00 );
  mdl_cjQtu1Re = slha.get_block_entry( "smeft", 82, 0.000000e+00 );
  mdl_cQb8 = slha.get_block_entry( "smeft", 81, 0.000000e+00 );
  mdl_cQb1 = slha.get_block_entry( "smeft", 80, 0.000000e+00 );
  mdl_cbj8 = slha.get_block_entry( "smeft", 79, 0.000000e+00 );
  mdl_cbj1 = slha.get_block_entry( "smeft", 78, 0.000000e+00 );
  mdl_cQd8 = slha.get_block_entry( "smeft", 77, 0.000000e+00 );
  mdl_cQd1 = slha.get_block_entry( "smeft", 76, 0.000000e+00 );
  mdl_cjd8 = slha.get_block_entry( "smeft", 75, 0.000000e+00 );
  mdl_cjd1 = slha.get_block_entry( "smeft", 74, 0.000000e+00 );
  mdl_cQt8 = slha.get_block_entry( "smeft", 73, 0.000000e+00 );
  mdl_cQt1 = slha.get_block_entry( "smeft", 72, 0.000000e+00 );
  mdl_ctj8 = slha.get_block_entry( "smeft", 71, 0.000000e+00 );
  mdl_ctj1 = slha.get_block_entry( "smeft", 70, 0.000000e+00 );
  mdl_cQu8 = slha.get_block_entry( "smeft", 69, 0.000000e+00 );
  mdl_cju8 = slha.get_block_entry( "smeft", 68, 0.000000e+00 );
  mdl_cQu1 = slha.get_block_entry( "smeft", 67, 0.000000e+00 );
  mdl_cju1 = slha.get_block_entry( "smeft", 66, 0.000000e+00 );
  mdl_cutbd8Re = slha.get_block_entry( "smeft", 65, 0.000000e+00 );
  mdl_cutbd1Re = slha.get_block_entry( "smeft", 64, 0.000000e+00 );
  mdl_cbu8 = slha.get_block_entry( "smeft", 63, 0.000000e+00 );
  mdl_ctd8 = slha.get_block_entry( "smeft", 62, 0.000000e+00 );
  mdl_ctb8 = slha.get_block_entry( "smeft", 61, 0.000000e+00 );
  mdl_cud8 = slha.get_block_entry( "smeft", 60, 0.000000e+00 );
  mdl_cbu1 = slha.get_block_entry( "smeft", 59, 0.000000e+00 );
  mdl_ctd1 = slha.get_block_entry( "smeft", 58, 0.000000e+00 );
  mdl_ctb1 = slha.get_block_entry( "smeft", 57, 0.000000e+00 );
  mdl_cud1 = slha.get_block_entry( "smeft", 56, 0.000000e+00 );
  mdl_cbd8 = slha.get_block_entry( "smeft", 55, 0.000000e+00 );
  mdl_cbd1 = slha.get_block_entry( "smeft", 54, 0.000000e+00 );
  mdl_cbb = slha.get_block_entry( "smeft", 53, 0.000000e+00 );
  mdl_cdd8 = slha.get_block_entry( "smeft", 52, 0.000000e+00 );
  mdl_cdd1 = slha.get_block_entry( "smeft", 51, 0.000000e+00 );
  mdl_ctu8 = slha.get_block_entry( "smeft", 50, 0.000000e+00 );
  mdl_ctu1 = slha.get_block_entry( "smeft", 49, 0.000000e+00 );
  mdl_ctt = slha.get_block_entry( "smeft", 48, 0.000000e+00 );
  mdl_cuu8 = slha.get_block_entry( "smeft", 47, 0.000000e+00 );
  mdl_cuu1 = slha.get_block_entry( "smeft", 46, 0.000000e+00 );
  mdl_cQQ8 = slha.get_block_entry( "smeft", 45, 0.000000e+00 );
  mdl_cQQ1 = slha.get_block_entry( "smeft", 44, 0.000000e+00 );
  mdl_cQj38 = slha.get_block_entry( "smeft", 43, 0.000000e+00 );
  mdl_cQj31 = slha.get_block_entry( "smeft", 42, 0.000000e+00 );
  mdl_cQj18 = slha.get_block_entry( "smeft", 41, 0.000000e+00 );
  mdl_cQj11 = slha.get_block_entry( "smeft", 40, 0.000000e+00 );
  mdl_cjj38 = slha.get_block_entry( "smeft", 39, 0.000000e+00 );
  mdl_cjj31 = slha.get_block_entry( "smeft", 38, 0.000000e+00 );
  mdl_cjj18 = slha.get_block_entry( "smeft", 37, 0.000000e+00 );
  mdl_cjj11 = slha.get_block_entry( "smeft", 36, 0.000000e+00 );
  mdl_cHtbRe = slha.get_block_entry( "smeft", 35, 0.000000e+00 );
  mdl_cHudRe = slha.get_block_entry( "smeft", 34, 0.000000e+00 );
  mdl_cHbq = slha.get_block_entry( "smeft", 33, 0.000000e+00 );
  mdl_cHd = slha.get_block_entry( "smeft", 32, 0.000000e+00 );
  mdl_cHt = slha.get_block_entry( "smeft", 31, 0.000000e+00 );
  mdl_cHu = slha.get_block_entry( "smeft", 30, 0.000000e+00 );
  mdl_cHQ3 = slha.get_block_entry( "smeft", 29, 0.000000e+00 );
  mdl_cHj3 = slha.get_block_entry( "smeft", 28, 0.000000e+00 );
  mdl_cHQ1 = slha.get_block_entry( "smeft", 27, 0.000000e+00 );
  mdl_cHj1 = slha.get_block_entry( "smeft", 26, 0.000000e+00 );
  mdl_cbBRe = slha.get_block_entry( "smeft", 25, 0.000000e+00 );
  mdl_cdBRe = slha.get_block_entry( "smeft", 24, 0.000000e+00 );
  mdl_cbWRe = slha.get_block_entry( "smeft", 23, 0.000000e+00 );
  mdl_cdWRe = slha.get_block_entry( "smeft", 22, 0.000000e+00 );
  mdl_cbGRe = slha.get_block_entry( "smeft", 21, 0.000000e+00 );
  mdl_cdGRe = slha.get_block_entry( "smeft", 20, 0.000000e+00 );
  mdl_ctBRe = slha.get_block_entry( "smeft", 19, 0.000000e+00 );
  mdl_cuBRe = slha.get_block_entry( "smeft", 18, 0.000000e+00 );
  mdl_ctWRe = slha.get_block_entry( "smeft", 17, 0.000000e+00 );
  mdl_cuWRe = slha.get_block_entry( "smeft", 16, 0.000000e+00 );
  mdl_ctGRe = slha.get_block_entry( "smeft", 15, 0.000000e+00 );
  mdl_cuGRe = slha.get_block_entry( "smeft", 14, 0.000000e+00 );
  mdl_cbHRe = slha.get_block_entry( "smeft", 13, 0.000000e+00 );
  mdl_cdHRe = slha.get_block_entry( "smeft", 12, 0.000000e+00 );
  mdl_ctHRe = slha.get_block_entry( "smeft", 11, 0.000000e+00 );
  mdl_cuHRe = slha.get_block_entry( "smeft", 10, 0.000000e+00 );
  mdl_cHWB = slha.get_block_entry( "smeft", 9, 0.000000e+00 );
  mdl_cHB = slha.get_block_entry( "smeft", 8, 0.000000e+00 );
  mdl_cHW = slha.get_block_entry( "smeft", 7, 0.000000e+00 );
  mdl_cHG = slha.get_block_entry( "smeft", 6, 0.000000e+00 );
  mdl_cHDD = slha.get_block_entry( "smeft", 5, 0.000000e+00 );
  mdl_cHbox = slha.get_block_entry( "smeft", 4, 0.000000e+00 );
  mdl_cH = slha.get_block_entry( "smeft", 3, 0.000000e+00 );
  mdl_cW = slha.get_block_entry( "smeft", 2, 0.000000e+00 );
  mdl_cG = slha.get_block_entry( "smeft", 1, 0.000000e+00 );
  mdl_MH = slha.get_block_entry( "mass", 25, 1.250900e+02 );
  mdl_MZ = slha.get_block_entry( "mass", 23, 9.118760e+01 );
  mdl_MTA = slha.get_block_entry( "mass", 15, 1.777000e+00 );
  mdl_MMU = slha.get_block_entry( "mass", 13, 1.056600e-01 );
  mdl_Me = slha.get_block_entry( "mass", 11, 5.110000e-04 );
  mdl_MT = slha.get_block_entry( "mass", 6, 1.727600e+02 );
  mdl_MB = slha.get_block_entry( "mass", 5, 4.180000e+00 );
  mdl_MC = slha.get_block_entry( "mass", 4, 1.270000e+00 );
  mdl_MS = slha.get_block_entry( "mass", 3, 9.300000e-02 );
  mdl_MU = slha.get_block_entry( "mass", 2, 2.160000e-03 );
  mdl_MD = slha.get_block_entry( "mass", 1, 4.670000e-03 );
  mdl_complexi = cxsmpl<double>( 0., 1. );
  mdl_cuH = mdl_cuHRe + mdl_cuHIm * mdl_complexi;
  mdl_ctHH = mdl_ctHRe + mdl_ctHIm * mdl_complexi;
  mdl_cdH = mdl_cdHRe + mdl_cdHIm * mdl_complexi;
  mdl_cbH = mdl_cbHRe + mdl_cbHIm * mdl_complexi;
  mdl_cuG = mdl_cuGRe + mdl_cuGIm * mdl_complexi;
  mdl_ctG = mdl_ctGRe + mdl_ctGIm * mdl_complexi;
  mdl_cuW = mdl_cuWRe + mdl_cuWIm * mdl_complexi;
  mdl_ctW = mdl_ctWRe + mdl_ctWIm * mdl_complexi;
  mdl_cuB = mdl_cuBRe + mdl_cuBIm * mdl_complexi;
  mdl_ctB = mdl_ctBRe + mdl_ctBIm * mdl_complexi;
  mdl_cdG = mdl_cdGRe + mdl_cdGIm * mdl_complexi;
  mdl_cbG = mdl_cbGRe + mdl_cbGIm * mdl_complexi;
  mdl_cdW = mdl_cdWRe + mdl_cdWIm * mdl_complexi;
  mdl_cbW = mdl_cbWRe + mdl_cbWIm * mdl_complexi;
  mdl_cdB = mdl_cdBRe + mdl_cdBIm * mdl_complexi;
  mdl_cbBB = mdl_cbBRe + mdl_cbBIm * mdl_complexi;
  mdl_cHud = mdl_cHudRe + mdl_cHudIm * mdl_complexi;
  mdl_cHtb = mdl_cHtbRe + mdl_cHtbIm * mdl_complexi;
  mdl_cutbd1 = mdl_cutbd1Re + mdl_cutbd1Im * mdl_complexi;
  mdl_cutbd8 = mdl_cutbd8Re + mdl_cutbd8Im * mdl_complexi;
  mdl_cjQtu1 = mdl_cjQtu1Re + mdl_cjQtu1Im * mdl_complexi;
  mdl_cjQtu8 = mdl_cjQtu8Re + mdl_cjQtu8Im * mdl_complexi;
  mdl_cjQbd1 = mdl_cjQbd1Re + mdl_cjQbd1Im * mdl_complexi;
  mdl_cjQbd8 = mdl_cjQbd8Re + mdl_cjQbd8Im * mdl_complexi;
  mdl_cjujd1 = mdl_cjujd1Re + mdl_cjujd1Im * mdl_complexi;
  mdl_cjujd8 = mdl_cjujd8Re + mdl_cjujd8Im * mdl_complexi;
  mdl_cjujd11 = mdl_cjujd11Re + mdl_cjujd11Im * mdl_complexi;
  mdl_cjujd81 = mdl_cjujd81Re + mdl_cjujd81Im * mdl_complexi;
  mdl_cQtjd1 = mdl_cQtjd1Re + mdl_cQtjd1Im * mdl_complexi;
  mdl_cQtjd8 = mdl_cQtjd8Re + mdl_cQtjd8Im * mdl_complexi;
  mdl_cjuQb1 = mdl_cjuQb1Re + mdl_cjuQb1Im * mdl_complexi;
  mdl_cjuQb8 = mdl_cjuQb8Re + mdl_cjuQb8Im * mdl_complexi;
  mdl_cQujb1 = mdl_cQujb1Re + mdl_cQujb1Im * mdl_complexi;
  mdl_cQujb8 = mdl_cQujb8Re + mdl_cQujb8Im * mdl_complexi;
  mdl_cjtQd1 = mdl_cjtQd1Re + mdl_cjtQd1Im * mdl_complexi;
  mdl_cjtQd8 = mdl_cjtQd8Re + mdl_cjtQd8Im * mdl_complexi;
  mdl_cQtQb1 = mdl_cQtQb1Re + mdl_cQtQb1Im * mdl_complexi;
  mdl_cQtQb8 = mdl_cQtQb8Re + mdl_cQtQb8Im * mdl_complexi;
  mdl_ceH = mdl_ceHRe + mdl_ceHIm * mdl_complexi;
  mdl_ceW = mdl_ceWRe + mdl_ceWIm * mdl_complexi;
  mdl_ceB = mdl_ceBRe + mdl_ceBIm * mdl_complexi;
  mdl_cledj = mdl_cledjRe + mdl_cledjIm * mdl_complexi;
  mdl_clebQ = mdl_clebQRe + mdl_clebQIm * mdl_complexi;
  mdl_cleju1 = mdl_cleju1Re + mdl_cleju1Im * mdl_complexi;
  mdl_cleju3 = mdl_cleju3Re + mdl_cleju3Im * mdl_complexi;
  mdl_cleQt1 = mdl_cleQt1Re + mdl_cleQt1Im * mdl_complexi;
  mdl_cleQt3 = mdl_cleQt3Re + mdl_cleQt3Im * mdl_complexi;
  mdl_MWsm = mdl_MW;
  mdl_MW__exp__2 = ( ( mdl_MW ) * ( mdl_MW ) );
  mdl_MZ__exp__2 = ( ( mdl_MZ ) * ( mdl_MZ ) );
  mdl_sqrt__2 = sqrt( 2. );
  mdl_nb__2__exp__0_25 = pow( 2., 0.25 );
  mdl_MH__exp__2 = ( ( mdl_MH ) * ( mdl_MH ) );
  mdl_sth2 = 1. - mdl_MW__exp__2 / mdl_MZ__exp__2;
  mdl_nb__10__exp___m_40 = pow( 10., -40. );
  mdl_propCorr = std::abs( mdl_linearPropCorrections ) / ( std::abs( mdl_linearPropCorrections ) + mdl_nb__10__exp___m_40 );
  mdl_MZ1 = mdl_MZ;
  mdl_MH1 = mdl_MH;
  mdl_MT1 = mdl_MT;
  mdl_WZ1 = mdl_WZ;
  mdl_WW1 = mdl_WW;
  mdl_WH1 = mdl_WH;
  mdl_WT1 = mdl_WT;
  mdl_cth = sqrt( 1. - mdl_sth2 );
  mdl_MW1 = mdl_MWsm;
  mdl_sqrt__sth2 = sqrt( mdl_sth2 );
  mdl_sth = mdl_sqrt__sth2;
  mdl_LambdaSMEFT__exp__2 = ( ( mdl_LambdaSMEFT ) * ( mdl_LambdaSMEFT ) );
  mdl_conjg__cbH = conj( mdl_cbH );
  mdl_conjg__ctHH = conj( mdl_ctHH );
  mdl_MT__exp__2 = ( ( mdl_MT ) * ( mdl_MT ) );
  mdl_MH__exp__6 = pow( mdl_MH, 6. );
  mdl_MWsm__exp__6 = pow( mdl_MWsm, 6. );
  mdl_MH__exp__4 = ( ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) * ( mdl_MH ) );
  mdl_MWsm__exp__4 = ( ( mdl_MWsm ) * ( mdl_MWsm ) * ( mdl_MWsm ) * ( mdl_MWsm ) );
  mdl_MWsm__exp__2 = ( ( mdl_MWsm ) * ( mdl_MWsm ) );
  mdl_MZ__exp__4 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
  mdl_MZ__exp__6 = pow( mdl_MZ, 6. );
  mdl_cth__exp__2 = ( ( mdl_cth ) * ( mdl_cth ) );
  mdl_sth__exp__2 = ( ( mdl_sth ) * ( mdl_sth ) );
  mdl_MB__exp__2 = ( ( mdl_MB ) * ( mdl_MB ) );
  mdl_MZ__exp__3 = ( ( mdl_MZ ) * ( mdl_MZ ) * ( mdl_MZ ) );
  mdl_sth__exp__4 = ( ( mdl_sth ) * ( mdl_sth ) * ( mdl_sth ) * ( mdl_sth ) );
  mdl_sth__exp__6 = pow( mdl_sth, 6. );
  mdl_sth__exp__3 = ( ( mdl_sth ) * ( mdl_sth ) * ( mdl_sth ) );
  mdl_sth__exp__5 = pow( mdl_sth, 5. );
  mdl_propCorr__exp__2 = ( ( mdl_propCorr ) * ( mdl_propCorr ) );
  mdl_propCorr__exp__3 = ( ( mdl_propCorr ) * ( mdl_propCorr ) * ( mdl_propCorr ) );
  mdl_propCorr__exp__4 = ( ( mdl_propCorr ) * ( mdl_propCorr ) * ( mdl_propCorr ) * ( mdl_propCorr ) );
  mdl_cth__exp__3 = ( ( mdl_cth ) * ( mdl_cth ) * ( mdl_cth ) );
  mdl_aEW = ( mdl_Gf * mdl_MW__exp__2 * ( 1. - mdl_MW__exp__2 / mdl_MZ__exp__2 ) * mdl_sqrt__2 ) / M_PI;
  mdl_sqrt__Gf = sqrt( mdl_Gf );
  mdl_vevhat = 1. / ( mdl_nb__2__exp__0_25 * mdl_sqrt__Gf );
  mdl_lam = ( mdl_Gf * mdl_MH__exp__2 ) / mdl_sqrt__2;
  mdl_sqrt__aEW = sqrt( mdl_aEW );
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt( M_PI );
  mdl_yb = ( mdl_ymb * mdl_sqrt__2 ) / mdl_vevhat;
  mdl_yc = ( mdl_ymc * mdl_sqrt__2 ) / mdl_vevhat;
  mdl_ydo = ( mdl_ymdo * mdl_sqrt__2 ) / mdl_vevhat;
  mdl_ye = ( mdl_yme * mdl_sqrt__2 ) / mdl_vevhat;
  mdl_ym = ( mdl_ymm * mdl_sqrt__2 ) / mdl_vevhat;
  mdl_ys = ( mdl_yms * mdl_sqrt__2 ) / mdl_vevhat;
  mdl_yt = ( mdl_ymt * mdl_sqrt__2 ) / mdl_vevhat;
  mdl_ytau = ( mdl_ymtau * mdl_sqrt__2 ) / mdl_vevhat;
  mdl_yup = ( mdl_ymup * mdl_sqrt__2 ) / mdl_vevhat;
  mdl_vevhat__exp__2 = ( ( mdl_vevhat ) * ( mdl_vevhat ) );
  mdl_dGf = ( ( 2. * mdl_cHl3 - mdl_cll1 ) * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2;
  mdl_dkH = ( ( mdl_cHbox - mdl_cHDD / 4. ) * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2;
  mdl_vevT = ( 1. + mdl_dGf / 2. ) * mdl_vevhat;
  mdl_g1 = mdl_ee / mdl_cth;
  mdl_gw = mdl_ee / mdl_sth;
  mdl_yb0 = ( 1. - mdl_dGf / 2. ) * mdl_yb + ( mdl_vevhat__exp__2 * mdl_conjg__cbH ) / ( 2. * mdl_LambdaSMEFT__exp__2 );
  mdl_yt0 = ( 1. - mdl_dGf / 2. ) * mdl_yt + ( mdl_vevhat__exp__2 * mdl_conjg__ctHH ) / ( 2. * mdl_LambdaSMEFT__exp__2 );
  mdl_ee__exp__2 = ( ( mdl_ee ) * ( mdl_ee ) );
  mdl_gHaa = ( mdl_ee__exp__2 * ( -1.75 + ( 4. * ( 0.3333333333333333 + ( 7. * mdl_MH__exp__2 ) / ( 360. * mdl_MT__exp__2 ) ) ) / 3. - ( 29. * mdl_MH__exp__6 ) / ( 16800. * mdl_MWsm__exp__6 ) - ( 19. * mdl_MH__exp__4 ) / ( 1680. * mdl_MWsm__exp__4 ) - ( 11. * mdl_MH__exp__2 ) / ( 120. * mdl_MWsm__exp__2 ) ) ) / ( 8. * ( ( M_PI ) * ( M_PI ) ) );
  mdl_gHza = ( mdl_ee__exp__2 * ( ( ( 0.4583333333333333 + ( 29. * mdl_MH__exp__6 ) / ( 100800. * mdl_MWsm__exp__6 ) + ( 19. * mdl_MH__exp__4 ) / ( 10080. * mdl_MWsm__exp__4 ) + ( 11. * mdl_MH__exp__2 ) / ( 720. * mdl_MWsm__exp__2 ) + ( mdl_MH__exp__4 * mdl_MZ__exp__2 ) / ( 2100. * mdl_MWsm__exp__6 ) + ( mdl_MH__exp__2 * mdl_MZ__exp__2 ) / ( 280. * mdl_MWsm__exp__4 ) + ( 7. * mdl_MZ__exp__2 ) / ( 180. * mdl_MWsm__exp__2 ) + ( 67. * mdl_MH__exp__2 * mdl_MZ__exp__4 ) / ( 100800. * mdl_MWsm__exp__6 ) + ( 53. * mdl_MZ__exp__4 ) / ( 10080. * mdl_MWsm__exp__4 ) + ( 43. * mdl_MZ__exp__6 ) / ( 50400. * mdl_MWsm__exp__6 ) - ( 31. * mdl_cth__exp__2 ) / ( 24. * mdl_sth__exp__2 ) - ( 29. * mdl_cth__exp__2 * mdl_MH__exp__6 ) / ( 20160. * mdl_MWsm__exp__6 * mdl_sth__exp__2 ) - ( 19. * mdl_cth__exp__2 * mdl_MH__exp__4 ) / ( 2016. * mdl_MWsm__exp__4 * mdl_sth__exp__2 ) - ( 11. * mdl_cth__exp__2 * mdl_MH__exp__2 ) / ( 144. * mdl_MWsm__exp__2 * mdl_sth__exp__2 ) - ( mdl_cth__exp__2 * mdl_MH__exp__4 * mdl_MZ__exp__2 ) / ( 560. * mdl_MWsm__exp__6 * mdl_sth__exp__2 ) - ( 31. * mdl_cth__exp__2 * mdl_MH__exp__2 * mdl_MZ__exp__2 ) / ( 2520. * mdl_MWsm__exp__4 * mdl_sth__exp__2 ) - ( mdl_cth__exp__2 * mdl_MZ__exp__2 ) / ( 9. * mdl_MWsm__exp__2 * mdl_sth__exp__2 ) - ( 43. * mdl_cth__exp__2 * mdl_MH__exp__2 * mdl_MZ__exp__4 ) / ( 20160. * mdl_MWsm__exp__6 * mdl_sth__exp__2 ) - ( 17. * mdl_cth__exp__2 * mdl_MZ__exp__4 ) / ( 1120. * mdl_MWsm__exp__4 * mdl_sth__exp__2 ) - ( 5. * mdl_cth__exp__2 * mdl_MZ__exp__6 ) / ( 2016. * mdl_MWsm__exp__6 * mdl_sth__exp__2 ) ) * mdl_sth ) / mdl_cth + ( ( 0.3333333333333333 + ( 7. * mdl_MH__exp__2 ) / ( 360. * mdl_MT__exp__2 ) + ( 11. * mdl_MZ__exp__2 ) / ( 360. * mdl_MT__exp__2 ) ) * ( 0.5 - ( 4. * mdl_sth__exp__2 ) / 3. ) ) / ( mdl_cth * mdl_sth ) ) ) / ( 4. * ( ( M_PI ) * ( M_PI ) ) );
  mdl_dMZ2 = ( ( mdl_cHDD / 2. + 2. * mdl_cHWB * mdl_cth * mdl_sth ) * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2;
  mdl_dMH2 = 2. * mdl_dkH - ( 3. * mdl_cH * mdl_vevhat__exp__2 ) / ( 2. * mdl_lam * mdl_LambdaSMEFT__exp__2 );
  mdl_dgw = -mdl_dGf / 2.;
  mdl_barlam = ( 1. - mdl_dGf - mdl_dMH2 ) * mdl_lam;
  mdl_dWT = 2. * mdl_WT * ( mdl_dgw + ( mdl_vevhat * ( mdl_ee * ( 3. * mdl_cHtbRe * mdl_MB * mdl_MT * mdl_MWsm__exp__2 + mdl_cHQ3 * ( ( ( mdl_MB__exp__2 - mdl_MT__exp__2 ) * ( mdl_MB__exp__2 - mdl_MT__exp__2 ) ) + ( mdl_MB__exp__2 + mdl_MT__exp__2 ) * mdl_MWsm__exp__2 - 2. * mdl_MWsm__exp__4 ) ) * mdl_vevhat + 6. * mdl_MWsm__exp__2 * ( mdl_ctWRe * mdl_MT * ( mdl_MB__exp__2 - mdl_MT__exp__2 + mdl_MWsm__exp__2 ) + mdl_cbWRe * mdl_MB * ( -mdl_MB__exp__2 + mdl_MT__exp__2 + mdl_MWsm__exp__2 ) ) * mdl_sth * mdl_sqrt__2 ) ) / ( mdl_ee * mdl_LambdaSMEFT__exp__2 * ( ( ( mdl_MB__exp__2 - mdl_MT__exp__2 ) * ( mdl_MB__exp__2 - mdl_MT__exp__2 ) ) + ( mdl_MB__exp__2 + mdl_MT__exp__2 ) * mdl_MWsm__exp__2 - 2. * mdl_MWsm__exp__4 ) ) );
  mdl_dWW = ( 2. * mdl_dgw + ( 2. * ( 2. * mdl_cHj3 + mdl_cHl3 ) * mdl_vevhat__exp__2 ) / ( 3. * mdl_LambdaSMEFT__exp__2 ) ) * mdl_WW;
  mdl_gwsh = ( mdl_ee * ( 1. + mdl_dgw - ( mdl_cHW * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 ) ) / mdl_sth;
  mdl_vev = ( 1. - ( 3. * mdl_cH * mdl_vevhat__exp__2 ) / ( 8. * mdl_lam * mdl_LambdaSMEFT__exp__2 ) ) * mdl_vevT;
  mdl_dg1 = ( -mdl_dGf - mdl_dMZ2 / mdl_sth__exp__2 ) / 2.;
  mdl_dWHc = mdl_yc / ( mdl_yc + mdl_nb__10__exp___m_40 ) * ( -0.02884 * mdl_dGf + ( ( 0.05768 * mdl_cHbox - 0.01442 * mdl_cHDD - 0.05768 * mdl_cuHRe ) * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 );
  mdl_dWHb = mdl_yb / ( mdl_yb + mdl_nb__10__exp___m_40 ) * ( mdl_vevhat__exp__2 * ( -1.1618 * mdl_cbHRe ) / ( mdl_LambdaSMEFT__exp__2 * ( mdl_yb + mdl_nb__10__exp___m_40 ) ) - 0.5809 * mdl_dGf + ( mdl_vevhat__exp__2 * ( 1.1618 * mdl_cHbox - 0.29045 * mdl_cHDD ) ) / ( mdl_LambdaSMEFT__exp__2 ) );
  mdl_dWHta = mdl_ytau / ( mdl_ytau + mdl_nb__10__exp___m_40 ) * ( -0.06256 * mdl_dGf + mdl_vevhat__exp__2 * ( -0.12512 * mdl_ceHRe + 0.12512 * mdl_cHbox - 0.03128 * mdl_cHDD ) / ( mdl_LambdaSMEFT__exp__2 ) );
  mdl_dWZ = mdl_WZ * ( -1. + ( 36. * mdl_cth * mdl_MB * mdl_MZ__exp__2 * mdl_sth * ( mdl_cbWRe * mdl_cth + mdl_cbBRe * mdl_sth ) * ( -3. + 4. * mdl_sth__exp__2 ) * mdl_vevhat * mdl_sqrt__2 * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_ee * mdl_LambdaSMEFT__exp__2 * ( 2. * mdl_MZ__exp__3 * ( 27. + 54. * mdl_dgw - 54. * ( 1. + mdl_dg1 + mdl_dgw ) * mdl_sth__exp__2 + 76. * ( 1. + 4. * mdl_dg1 - 2. * mdl_dgw ) * mdl_sth__exp__4 + 152. * ( -mdl_dg1 + mdl_dgw ) * mdl_sth__exp__6 ) + mdl_MZ__exp__2 * ( 9. + 18. * mdl_dgw - 6. * ( 2. + mdl_dg1 + 3. * mdl_dgw ) * mdl_sth__exp__2 + 8. * ( 1. + 4. * mdl_dg1 - 2. * mdl_dgw ) * mdl_sth__exp__4 + 16. * ( -mdl_dg1 + mdl_dgw ) * mdl_sth__exp__6 ) * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_MB__exp__2 * ( -9. - 18. * mdl_dgw - 6. * ( 4. + 11. * mdl_dg1 - 3. * mdl_dgw ) * mdl_sth__exp__2 + 16. * ( 1. + 4. * mdl_dg1 - 2. * mdl_dgw ) * mdl_sth__exp__4 + 32. * ( -mdl_dg1 + mdl_dgw ) * mdl_sth__exp__6 ) * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) ) + 2. * mdl_ee * mdl_vevhat__exp__2 * ( 36. * mdl_cHj3 * mdl_MZ__exp__3 + 18. * mdl_cHl3 * mdl_MZ__exp__3 + 9. * ( 3. * mdl_cHbq - mdl_cHQ1 - mdl_cHQ3 ) * mdl_MB__exp__2 * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + 9. * mdl_cHQ1 * mdl_MZ__exp__2 * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + 9. * mdl_cHQ3 * mdl_MZ__exp__2 * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + 3. * mdl_cHWB * mdl_cth * ( -7. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) * mdl_sth * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + 8. * mdl_cHWB * mdl_cth * mdl_sth__exp__3 * ( 2. * mdl_MB__exp__2 * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_MZ__exp__2 * ( 19. * mdl_MZ + sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) ) ) - 8. * mdl_cHWB * mdl_cth * mdl_sth__exp__5 * ( 2. * mdl_MB__exp__2 * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_MZ__exp__2 * ( 19. * mdl_MZ + sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) ) ) - 6. * mdl_sth__exp__2 * ( 2. * ( mdl_cHbq + mdl_cHQ1 + mdl_cHQ3 ) * mdl_MB__exp__2 * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_MZ__exp__2 * ( ( 2. * mdl_cHd + 3. * mdl_cHe - 2. * mdl_cHj1 + 3. * ( 2. * mdl_cHj3 + mdl_cHl1 + mdl_cHl3 ) - 4. * mdl_cHu ) * mdl_MZ + ( mdl_cHbq + mdl_cHQ1 + mdl_cHQ3 ) * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) ) ) ) ) / ( mdl_ee * mdl_LambdaSMEFT__exp__2 * ( 2. * mdl_MZ__exp__3 * ( 27. - 54. * mdl_sth__exp__2 + 76. * mdl_sth__exp__4 ) + mdl_MZ__exp__2 * ( 9. - 12. * mdl_sth__exp__2 + 8. * mdl_sth__exp__4 ) * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) + mdl_MB__exp__2 * ( -9. - 24. * mdl_sth__exp__2 + 16. * mdl_sth__exp__4 ) * sqrt( -4. * mdl_MB__exp__2 + mdl_MZ__exp__2 ) ) ) );
  mdl_g1sh = ( mdl_ee * ( 1. + mdl_dg1 - ( mdl_cHB * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 ) ) / mdl_cth;
  mdl_ee__exp__3 = ( ( mdl_ee ) * ( mdl_ee ) * ( mdl_ee ) );
  mdl_vevhat__exp__3 = ( ( mdl_vevhat ) * ( mdl_vevhat ) * ( mdl_vevhat ) );
  // BSM parameters that do not depend on alphaS but are needed in the computation of alphaS-dependent couplings;
  mdl_bsmIndepParam[0] = mdl_WH;
  mdl_bsmIndepParam[1] = mdl_LambdaSMEFT;
  mdl_bsmIndepParam[2] = mdl_cHl3;
  mdl_bsmIndepParam[3] = mdl_cHj3;
  mdl_bsmIndepParam[4] = mdl_cHWB;
  mdl_bsmIndepParam[5] = mdl_cHB;
  mdl_bsmIndepParam[6] = mdl_cHW;
  mdl_bsmIndepParam[7] = mdl_cHG;
  mdl_bsmIndepParam[8] = mdl_cH;
  mdl_bsmIndepParam[9] = mdl_MH;
  mdl_bsmIndepParam[10] = mdl_MT;
  mdl_bsmIndepParam[11] = mdl_MH__exp__2;
  mdl_bsmIndepParam[12] = mdl_cth;
  mdl_bsmIndepParam[13] = mdl_sth;
  mdl_bsmIndepParam[14] = mdl_LambdaSMEFT__exp__2;
  mdl_bsmIndepParam[15] = mdl_MT__exp__2;
  mdl_bsmIndepParam[16] = mdl_sth__exp__2;
  mdl_bsmIndepParam[17] = mdl_vevhat;
  mdl_bsmIndepParam[18] = mdl_vevhat__exp__2;
  mdl_bsmIndepParam[19] = mdl_dGf;
  mdl_bsmIndepParam[20] = mdl_dkH;
  mdl_bsmIndepParam[21] = mdl_gHaa;
  mdl_bsmIndepParam[22] = mdl_gHza;
  mdl_bsmIndepParam[23] = mdl_dgw;
  mdl_bsmIndepParam[24] = mdl_dWW;
  mdl_bsmIndepParam[25] = mdl_vev;
  mdl_bsmIndepParam[26] = mdl_dWHc;
  mdl_bsmIndepParam[27] = mdl_dWHb;
  mdl_bsmIndepParam[28] = mdl_dWHta;
}

void
Parameters_SMEFTsim_topU3l_MwScheme_UFO::setIndependentCouplings()
{
  // (none)
}

/*
void
Parameters_SMEFTsim_topU3l_MwScheme_UFO::setDependentParameters() // now computed event-by-event (running alphas #373)
{
  mdl_sqrt__aS = sqrt( aS );
  G = 2. * mdl_sqrt__aS * sqrt( M_PI );
  mdl_gHgg2 = ( -7. * aS ) / ( 720. * M_PI );
  mdl_gHgg4 = aS / ( 360. * M_PI );
  mdl_gHgg5 = aS / ( 20. * M_PI );
  mdl_G__exp__2 = ( ( G ) * ( G ) );
  mdl_gHgg1 = mdl_G__exp__2 / ( 48. * ( ( M_PI ) * ( M_PI ) ) );
  mdl_gHgg3 = ( aS * G ) / ( 60. * M_PI );
  mdl_G__exp__3 = ( ( G ) * ( G ) * ( G ) );
  mdl_dWH = mdl_WH * ( -0.24161 * mdl_dGf + 0.96644 * mdl_dgw + 0.4832199999999999 * mdl_dkH - 0.11186509426655467 * mdl_dWW + ( 0.36410378449238195 * mdl_cHj3 * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 + ( 0.17608307708657747 * mdl_cHl3 * mdl_vevhat__exp__2 ) / mdl_LambdaSMEFT__exp__2 + ( 0.1636 * mdl_cHG * mdl_MT__exp__2 * mdl_vevhat__exp__2 ) / ( mdl_LambdaSMEFT__exp__2 * ( -0.5 * mdl_gHgg2 * mdl_MH__exp__2 + mdl_gHgg1 * mdl_MT__exp__2 ) ) + ( mdl_cHW * ( -0.35937785117066967 * mdl_gHaa * mdl_gHza + 0.006164 * mdl_cth * mdl_gHaa * mdl_sth + 0.00454 * mdl_gHza * mdl_sth__exp__2 ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + ( mdl_cHWB * ( -0.00454 * mdl_cth * mdl_gHza * mdl_sth + mdl_gHaa * ( -0.0030819999999999997 + 0.006163999999999999 * mdl_sth__exp__2 ) ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + ( mdl_cHB * ( -0.006163999999999999 * mdl_cth * mdl_gHaa * mdl_sth - 0.00454 * mdl_gHza * ( -1. + mdl_sth__exp__2 ) ) * mdl_vevhat__exp__2 ) / ( mdl_gHaa * mdl_gHza * mdl_LambdaSMEFT__exp__2 ) + mdl_dWHc + mdl_dWHb + mdl_dWHta );
}

void
Parameters_SMEFTsim_topU3l_MwScheme_UFO::setDependentCouplings() // now computed event-by-event (running alphas #373)
{
  GC_7 = G;
  GC_6 = -( mdl_complexi * G );
  GC_8 = mdl_complexi * mdl_G__exp__2;
}
*/

#endif

// Routines for printing out parameters
void
Parameters_SMEFTsim_topU3l_MwScheme_UFO::printIndependentParameters()
{
  std::cout << "SMEFTsim_topU3l_MwScheme_UFO model parameters independent of event kinematics:" << std::endl;
  std::cout << "(Warning: aS in the runcard is ignored because event-by-event Gs are hardcoded or retrieved from Fortran)" << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WZ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WZ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WT = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WT << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymtau = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymtau << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yme = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yme << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymt = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymt << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymb = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymb << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymc = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymc << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yms = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yms << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymup = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymup << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ymdo = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ymdo << std::endl;
  std::cout << std::setw( 20 ) << "mdl_linearPropCorrections = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_linearPropCorrections << std::endl;
  //std::cout << std::setw( 20 ) << "aS = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << aS << std::endl; // now retrieved event-by-event (as G) from Fortran (running alphas #373)
  std::cout << std::setw( 20 ) << "mdl_Gf = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Gf << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_LambdaSMEFT = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_LambdaSMEFT << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleQt3Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleQt3Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleQt1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleQt1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleju3Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleju3Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleju1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleju1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_clebQIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_clebQIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cledjIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cledjIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ceBIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ceBIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ceWIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ceWIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ceHIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ceHIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtQb8Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtQb8Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtQb1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtQb1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjtQd8Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjtQd8Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjtQd1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjtQd1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQujb8Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQujb8Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQujb1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQujb1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjuQb8Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjuQb8Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjuQb1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjuQb1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtjd8Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtjd8Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtjd1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtjd1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd81Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd81Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd11Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd11Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd8Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd8Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQbd8Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQbd8Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQbd1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQbd1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQtu8Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQtu8Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQtu1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQtu1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cutbd8Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cutbd8Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cutbd1Im = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cutbd1Im << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHtbIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHtbIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHudIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHudIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbHIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbHIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdHIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdHIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctHIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctHIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuHIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuHIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbBIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbBIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdBIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdBIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbWIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbWIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdWIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdWIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbGIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbGIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdGIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdGIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctBIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctBIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuBIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuBIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctWIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctWIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuWIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuWIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctGIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctGIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuGIm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuGIm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHWBtil = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHWBtil << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHBtil = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHBtil << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHWtil = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHWtil << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHGtil = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHGtil << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cWtil = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cWtil << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cGtil = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cGtil << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleQt3Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleQt3Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleju3Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleju3Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleQt1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleQt1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleju1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleju1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_clebQRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_clebQRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cledjRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cledjRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cle = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cle << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbl = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbl << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cld = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cld << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctl = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctl << std::endl;
  std::cout << std::setw( 20 ) << "mdl_clu = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_clu << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cje = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cje << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ced = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ced << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cte = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cte << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ceu = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ceu << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cee = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cee << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQl3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQl3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQl1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQl1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_clj3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_clj3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_clj1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_clj1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cll1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cll1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cll = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cll << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHl3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHl3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHl1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHl1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ceBRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ceBRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ceWRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ceWRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ceHRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ceHRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtQb8Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtQb8Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtQb1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtQb1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjtQd8Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjtQd8Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjtQd1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjtQd1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQujb8Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQujb8Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQujb1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQujb1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjuQb8Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjuQb8Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjuQb1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjuQb1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtjd8Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtjd8Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtjd1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtjd1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd81Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd81Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd11Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd11Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd8Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd8Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQbd8Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQbd8Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQbd1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQbd1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQtu8Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQtu8Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQtu1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQtu1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQb8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQb8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQb1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQb1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbj8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbj8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbj1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbj1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQd8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQd8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQd1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQd1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjd8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjd8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjd1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjd1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQt8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQt8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQt1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQt1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctj8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctj8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctj1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctj1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQu8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQu8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cju8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cju8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQu1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQu1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cju1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cju1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cutbd8Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cutbd8Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cutbd1Re = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cutbd1Re << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbu8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbu8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctd8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctd8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctb8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctb8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cud8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cud8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbu1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbu1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctd1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctd1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctb1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctb1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cud1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cud1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbd8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbd8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbd1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbd1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbb = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbb << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdd8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdd8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdd1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdd1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctu8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctu8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctu1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctu1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctt = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctt << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuu8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuu8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuu1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuu1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQQ8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQQ8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQQ1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQQ1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQj38 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQj38 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQj31 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQj31 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQj18 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQj18 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQj11 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQj11 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjj38 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjj38 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjj31 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjj31 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjj18 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjj18 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjj11 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjj11 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHtbRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHtbRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHudRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHudRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHbq = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHbq << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHd = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHd << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHt = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHt << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHu = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHu << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHQ3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHQ3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHj3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHj3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHQ1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHQ1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHj1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHj1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbBRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbBRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdBRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdBRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbWRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbWRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdWRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdWRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbGRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbGRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdGRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdGRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctBRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctBRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuBRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuBRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctWRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctWRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuWRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuWRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctGRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctGRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuGRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuGRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbHRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbHRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdHRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdHRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctHRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctHRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuHRe = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuHRe << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHWB = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHWB << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHB = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHB << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHG = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHG << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHDD = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHDD << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHbox = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHbox << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cG = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cG << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MTA = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MTA << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MMU = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MMU << std::endl;
  std::cout << std::setw( 20 ) << "mdl_Me = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_Me << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MB = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MB << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MC = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MC << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MS = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MS << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MU = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MU << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MD = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MD << std::endl;
  std::cout << std::setw( 20 ) << "mdl_complexi = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_complexi << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctHH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctHH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuG = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuG << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctG = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctG << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cuB = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cuB << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ctB = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ctB << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdG = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdG << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbG = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbG << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cdB = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cdB << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cbBB = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cbBB << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHud = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHud << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cHtb = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cHtb << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cutbd1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cutbd1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cutbd8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cutbd8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQtu1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQtu1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQtu8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQtu8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQbd1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQbd1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjQbd8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjQbd8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd11 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd11 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjujd81 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjujd81 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtjd1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtjd1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtjd8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtjd8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjuQb1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjuQb1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjuQb8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjuQb8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQujb1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQujb1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQujb8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQujb8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjtQd1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjtQd1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cjtQd8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cjtQd8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtQb1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtQb1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cQtQb8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cQtQb8 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ceH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ceH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ceW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ceW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ceB = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ceB << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cledj = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cledj << std::endl;
  std::cout << std::setw( 20 ) << "mdl_clebQ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_clebQ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleju1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleju1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleju3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleju3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleQt1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleQt1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cleQt3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cleQt3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MWsm = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MWsm << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_nb__2__exp__0_25 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_nb__2__exp__0_25 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sth2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sth2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_nb__10__exp___m_40 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_nb__10__exp___m_40 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_propCorr = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_propCorr << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WZ1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WZ1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WW1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WW1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WH1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WH1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_WT1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_WT1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cth = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cth << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MW1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MW1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__sth2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__sth2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sth = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sth << std::endl;
  std::cout << std::setw( 20 ) << "mdl_LambdaSMEFT__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_LambdaSMEFT__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_conjg__cbH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_conjg__cbH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_conjg__ctHH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_conjg__ctHH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MT__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MT__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MWsm__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MWsm__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MH__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MH__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MWsm__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MWsm__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MWsm__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MWsm__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cth__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cth__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sth__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sth__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MB__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MB__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_MZ__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_MZ__exp__3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sth__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sth__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sth__exp__6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sth__exp__6 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sth__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sth__exp__3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sth__exp__5 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sth__exp__5 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_propCorr__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_propCorr__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_propCorr__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_propCorr__exp__3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_propCorr__exp__4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_propCorr__exp__4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_cth__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_cth__exp__3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_aEW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_aEW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__Gf = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__Gf << std::endl;
  std::cout << std::setw( 20 ) << "mdl_vevhat = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_vevhat << std::endl;
  std::cout << std::setw( 20 ) << "mdl_lam = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_lam << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__aEW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__aEW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yb = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yb << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yc = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yc << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ydo = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ydo << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ye = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ye << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ym = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ym << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ys = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ys << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yt = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yt << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ytau = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ytau << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yup = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yup << std::endl;
  std::cout << std::setw( 20 ) << "mdl_vevhat__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_vevhat__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dGf = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dGf << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dkH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dkH << std::endl;
  std::cout << std::setw( 20 ) << "mdl_vevT = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_vevT << std::endl;
  std::cout << std::setw( 20 ) << "mdl_g1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_g1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gw = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gw << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yb0 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yb0 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_yt0 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_yt0 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gHaa = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gHaa << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gHza = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gHza << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dMZ2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dMZ2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dMH2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dMH2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dgw = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dgw << std::endl;
  std::cout << std::setw( 20 ) << "mdl_barlam = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_barlam << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dWT = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dWT << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dWW = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dWW << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gwsh = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gwsh << std::endl;
  std::cout << std::setw( 20 ) << "mdl_vev = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_vev << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dg1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dg1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dWHc = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dWHc << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dWHb = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dWHb << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dWHta = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dWHta << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dWZ = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dWZ << std::endl;
  std::cout << std::setw( 20 ) << "mdl_g1sh = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_g1sh << std::endl;
  std::cout << std::setw( 20 ) << "mdl_ee__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_ee__exp__3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_vevhat__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_vevhat__exp__3 << std::endl;
}

void
Parameters_SMEFTsim_topU3l_MwScheme_UFO::printIndependentCouplings()
{
  std::cout << "SMEFTsim_topU3l_MwScheme_UFO model couplings independent of event kinematics:" << std::endl;
  // (none)
}

/*
void
Parameters_SMEFTsim_topU3l_MwScheme_UFO::printDependentParameters() // now computed event-by-event (running alphas #373)
{
  std::cout << "SMEFTsim_topU3l_MwScheme_UFO model parameters dependent on event kinematics:" << std::endl;
  std::cout << std::setw( 20 ) << "mdl_sqrt__aS = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_sqrt__aS << std::endl;
  std::cout << std::setw( 20 ) << "G = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << G << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gHgg2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gHgg2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gHgg4 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gHgg4 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gHgg5 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gHgg5 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G__exp__2 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G__exp__2 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gHgg1 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gHgg1 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_gHgg3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_gHgg3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_G__exp__3 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_G__exp__3 << std::endl;
  std::cout << std::setw( 20 ) << "mdl_dWH = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << mdl_dWH << std::endl;
}

void
Parameters_SMEFTsim_topU3l_MwScheme_UFO::printDependentCouplings() // now computed event-by-event (running alphas #373)
{
  std::cout << "SMEFTsim_topU3l_MwScheme_UFO model couplings dependent on event kinematics:" << std::endl;
  std::cout << std::setw( 20 ) << "GC_7 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_7 << std::endl;
  std::cout << std::setw( 20 ) << "GC_6 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_6 << std::endl;
  std::cout << std::setw( 20 ) << "GC_8 = " << std::setiosflags( std::ios::scientific ) << std::setw( 10 ) << GC_8 << std::endl;
}
*/
