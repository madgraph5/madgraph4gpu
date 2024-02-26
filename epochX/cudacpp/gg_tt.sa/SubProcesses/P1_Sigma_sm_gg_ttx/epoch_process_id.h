// Copyright (C) 2020-2023 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Oct 2021) for the MG5aMC CUDACPP plugin.
// Further modified by: A. Valassi (2021-2023) for the MG5aMC CUDACPP plugin.

#ifndef EPOCH_PROCESS_ID_H
#define EPOCH_PROCESS_ID_H 1

// No need to indicate EPOCHX_ any longer for auto-generated code
// However, keep the name of the file as it may be useful again for new manual developments
#define MG_EPOCH_PROCESS_ID SIGMA_SM_GG_TTX

// For simplicity, define here the name of the process-dependent reference file for tests
#define MG_EPOCH_REFERENCE_FILE_NAME "../../test/ref/dump_CPUTest.Sigma_sm_gg_ttx.txt"

#endif // EPOCH_PROCESS_ID_H
