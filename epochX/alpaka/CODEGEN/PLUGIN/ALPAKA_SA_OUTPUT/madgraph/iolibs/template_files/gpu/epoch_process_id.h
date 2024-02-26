#ifndef EPOCH_PROCESS_ID_H
#define EPOCH_PROCESS_ID_H 1

// No need to indicate EPOCHX_ any longer for auto-generated code
// However, keep the name of the file as it may be useful again for new manual developments
#define MG_EPOCH_PROCESS_ID %(processid_uppercase)s

// For simplicity, define here the name of the process-dependent reference file for tests
#define MG_EPOCH_REFERENCE_FILE_NAME "../../../../../test/ref/dump_CPUTest.%(processid)s.txt"

#endif // EPOCH_PROCESS_ID_H
