// Copyright (C) 2020-2024 CERN and UCLouvain.
// Licensed under the GNU Lesser General Public License (version 3 or later).
// Created by: A. Valassi (Dec 2022) for the MG5aMC CUDACPP plugin.
// Further modified by: O. Mattelaer, A. Valassi (2022-2024) for the MG5aMC CUDACPP plugin.

#ifndef COLORAMPS_H
#define COLORAMPS_H 1

// Note: strictly speaking the check '#ifdef MGONGPU_SUPPORTS_MULTICHANNEL' is not needed here,
// because coloramps.h is not included otherwise, but adding it does not harm and makes the code clearer

#ifdef MGONGPU_SUPPORTS_MULTICHANNEL /* clang-format off */

namespace mgOnGpu
{
  // Summary of numbering and indexing conventions for the relevant concepts (see issue #826 and PR #852)
  // - Diagram number (no variable) in [1, N_diagrams]: all values are allowed (N_diagrams distinct values)
  //   => this number is displayed for information before each block of code in CPPProcess.cc
  // - Channel number ("channelId" in C, CHANNEL_ID in F) in [1, N_diagrams]: not all values are allowed (N_config <= N_diagrams distinct values)
  //   => this number (with F indexing as in ps/pdf output) is passed around as an API argument between cudacpp functions
  //   Note: the old API passes around a single CHANNEL_ID (and uses CHANNEL_ID=0 to indicate no-multichannel mode, but this is not used in coloramps.h),
  //   while the new API passes around an array of CHANNEL_ID's (and uses a NULL array pointer to indicate no-multichannel mode)
  // - Channel number in C indexing: "channelID - 1"
  //   => this number (with C indexing) is used as the index of the channel2iconfig array below
  // - Config number ("iconfig" in C, ICONFIG in F) in [1, N_config]: all values are allowed (N_config <= N_diagrams distinct values)
  // - Config number in C indexing: "iconfig - 1"
  //   => this number (with C indexing) is used as the index of the icolamp array below

  // Map channel to iconfig (e.g. "iconfig = channel2iconfig[channelId - 1]": input index uses C indexing, output index uses F indexing)
  // Note: iconfig=-1 indicates channels/diagrams with no associated iconfig for single-diagram enhancement in the MadEvent sampling algorithm (presence of 4-point interaction?)
  // This array has N_diagrams elements, but only N_config <= N_diagrams valid values (iconfig>0)
  __device__ constexpr int channel2iconfig[48] = { // note: a trailing comma in the initializer list is allowed
      1, // CHANNEL_ID=1   i.e. DIAGRAM=1   --> ICONFIG=1
      2, // CHANNEL_ID=2   i.e. DIAGRAM=2   --> ICONFIG=2
      3, // CHANNEL_ID=3   i.e. DIAGRAM=3   --> ICONFIG=3
      4, // CHANNEL_ID=4   i.e. DIAGRAM=4   --> ICONFIG=4
      5, // CHANNEL_ID=5   i.e. DIAGRAM=5   --> ICONFIG=5
      6, // CHANNEL_ID=6   i.e. DIAGRAM=6   --> ICONFIG=6
      7, // CHANNEL_ID=7   i.e. DIAGRAM=7   --> ICONFIG=7
      8, // CHANNEL_ID=8   i.e. DIAGRAM=8   --> ICONFIG=8
      9, // CHANNEL_ID=9   i.e. DIAGRAM=9   --> ICONFIG=9
     10, // CHANNEL_ID=10  i.e. DIAGRAM=10  --> ICONFIG=10
     11, // CHANNEL_ID=11  i.e. DIAGRAM=11  --> ICONFIG=11
     12, // CHANNEL_ID=12  i.e. DIAGRAM=12  --> ICONFIG=12
     13, // CHANNEL_ID=13  i.e. DIAGRAM=13  --> ICONFIG=13
     14, // CHANNEL_ID=14  i.e. DIAGRAM=14  --> ICONFIG=14
     15, // CHANNEL_ID=15  i.e. DIAGRAM=15  --> ICONFIG=15
     16, // CHANNEL_ID=16  i.e. DIAGRAM=16  --> ICONFIG=16
     17, // CHANNEL_ID=17  i.e. DIAGRAM=17  --> ICONFIG=17
     18, // CHANNEL_ID=18  i.e. DIAGRAM=18  --> ICONFIG=18
     19, // CHANNEL_ID=19  i.e. DIAGRAM=19  --> ICONFIG=19
     20, // CHANNEL_ID=20  i.e. DIAGRAM=20  --> ICONFIG=20
     21, // CHANNEL_ID=21  i.e. DIAGRAM=21  --> ICONFIG=21
     22, // CHANNEL_ID=22  i.e. DIAGRAM=22  --> ICONFIG=22
     23, // CHANNEL_ID=23  i.e. DIAGRAM=23  --> ICONFIG=23
     24, // CHANNEL_ID=24  i.e. DIAGRAM=24  --> ICONFIG=24
     25, // CHANNEL_ID=25  i.e. DIAGRAM=25  --> ICONFIG=25
     26, // CHANNEL_ID=26  i.e. DIAGRAM=26  --> ICONFIG=26
     27, // CHANNEL_ID=27  i.e. DIAGRAM=27  --> ICONFIG=27
     28, // CHANNEL_ID=28  i.e. DIAGRAM=28  --> ICONFIG=28
     29, // CHANNEL_ID=29  i.e. DIAGRAM=29  --> ICONFIG=29
     30, // CHANNEL_ID=30  i.e. DIAGRAM=30  --> ICONFIG=30
     31, // CHANNEL_ID=31  i.e. DIAGRAM=31  --> ICONFIG=31
     32, // CHANNEL_ID=32  i.e. DIAGRAM=32  --> ICONFIG=32
     33, // CHANNEL_ID=33  i.e. DIAGRAM=33  --> ICONFIG=33
     34, // CHANNEL_ID=34  i.e. DIAGRAM=34  --> ICONFIG=34
     35, // CHANNEL_ID=35  i.e. DIAGRAM=35  --> ICONFIG=35
     36, // CHANNEL_ID=36  i.e. DIAGRAM=36  --> ICONFIG=36
     37, // CHANNEL_ID=37  i.e. DIAGRAM=37  --> ICONFIG=37
     38, // CHANNEL_ID=38  i.e. DIAGRAM=38  --> ICONFIG=38
     39, // CHANNEL_ID=39  i.e. DIAGRAM=39  --> ICONFIG=39
     40, // CHANNEL_ID=40  i.e. DIAGRAM=40  --> ICONFIG=40
     41, // CHANNEL_ID=41  i.e. DIAGRAM=41  --> ICONFIG=41
     42, // CHANNEL_ID=42  i.e. DIAGRAM=42  --> ICONFIG=42
     43, // CHANNEL_ID=43  i.e. DIAGRAM=43  --> ICONFIG=43
     44, // CHANNEL_ID=44  i.e. DIAGRAM=44  --> ICONFIG=44
     45, // CHANNEL_ID=45  i.e. DIAGRAM=45  --> ICONFIG=45
     46, // CHANNEL_ID=46  i.e. DIAGRAM=46  --> ICONFIG=46
     47, // CHANNEL_ID=47  i.e. DIAGRAM=47  --> ICONFIG=47
     48, // CHANNEL_ID=48  i.e. DIAGRAM=48  --> ICONFIG=48
  };

  // Map iconfig to the mask of allowed colors (e.g. "colormask = icolamp[iconfig - 1]": input index uses C indexing)
  // This array has N_config <= N_diagrams elements
  __device__ constexpr bool icolamp[48][4] = { // note: a trailing comma in the initializer list is allowed
    {  true, false, false, false }, // ICONFIG=1   <-- CHANNEL_ID=1
    {  true, false, false, false }, // ICONFIG=2   <-- CHANNEL_ID=2
    {  true, false, false, false }, // ICONFIG=3   <-- CHANNEL_ID=3
    {  true, false, false, false }, // ICONFIG=4   <-- CHANNEL_ID=4
    {  true, false, false, false }, // ICONFIG=5   <-- CHANNEL_ID=5
    {  true, false, false, false }, // ICONFIG=6   <-- CHANNEL_ID=6
    {  true, false, false, false }, // ICONFIG=7   <-- CHANNEL_ID=7
    {  true, false, false, false }, // ICONFIG=8   <-- CHANNEL_ID=8
    { false, false, false,  true }, // ICONFIG=9   <-- CHANNEL_ID=9
    { false, false, false,  true }, // ICONFIG=10  <-- CHANNEL_ID=10
    { false, false, false,  true }, // ICONFIG=11  <-- CHANNEL_ID=11
    { false, false, false,  true }, // ICONFIG=12  <-- CHANNEL_ID=12
    { false, false, false,  true }, // ICONFIG=13  <-- CHANNEL_ID=13
    { false, false, false,  true }, // ICONFIG=14  <-- CHANNEL_ID=14
    { false, false, false,  true }, // ICONFIG=15  <-- CHANNEL_ID=15
    { false, false, false,  true }, // ICONFIG=16  <-- CHANNEL_ID=16
    {  true, false, false, false }, // ICONFIG=17  <-- CHANNEL_ID=17
    {  true, false, false, false }, // ICONFIG=18  <-- CHANNEL_ID=18
    {  true, false, false, false }, // ICONFIG=19  <-- CHANNEL_ID=19
    {  true, false, false, false }, // ICONFIG=20  <-- CHANNEL_ID=20
    {  true, false, false, false }, // ICONFIG=21  <-- CHANNEL_ID=21
    {  true, false, false, false }, // ICONFIG=22  <-- CHANNEL_ID=22
    {  true, false, false, false }, // ICONFIG=23  <-- CHANNEL_ID=23
    {  true, false, false, false }, // ICONFIG=24  <-- CHANNEL_ID=24
    { false, false, false,  true }, // ICONFIG=25  <-- CHANNEL_ID=25
    { false, false, false,  true }, // ICONFIG=26  <-- CHANNEL_ID=26
    { false, false, false,  true }, // ICONFIG=27  <-- CHANNEL_ID=27
    { false, false, false,  true }, // ICONFIG=28  <-- CHANNEL_ID=28
    { false, false, false,  true }, // ICONFIG=29  <-- CHANNEL_ID=29
    { false, false, false,  true }, // ICONFIG=30  <-- CHANNEL_ID=30
    { false, false, false,  true }, // ICONFIG=31  <-- CHANNEL_ID=31
    { false, false, false,  true }, // ICONFIG=32  <-- CHANNEL_ID=32
    {  true, false, false,  true }, // ICONFIG=33  <-- CHANNEL_ID=33
    {  true, false, false,  true }, // ICONFIG=34  <-- CHANNEL_ID=34
    { false, false, false,  true }, // ICONFIG=35  <-- CHANNEL_ID=35
    {  true, false, false, false }, // ICONFIG=36  <-- CHANNEL_ID=36
    {  true, false, false,  true }, // ICONFIG=37  <-- CHANNEL_ID=37
    {  true, false, false,  true }, // ICONFIG=38  <-- CHANNEL_ID=38
    { false, false, false,  true }, // ICONFIG=39  <-- CHANNEL_ID=39
    {  true, false, false, false }, // ICONFIG=40  <-- CHANNEL_ID=40
    {  true, false, false,  true }, // ICONFIG=41  <-- CHANNEL_ID=41
    {  true, false, false,  true }, // ICONFIG=42  <-- CHANNEL_ID=42
    {  true, false, false, false }, // ICONFIG=43  <-- CHANNEL_ID=43
    { false, false, false,  true }, // ICONFIG=44  <-- CHANNEL_ID=44
    {  true, false, false,  true }, // ICONFIG=45  <-- CHANNEL_ID=45
    {  true, false, false,  true }, // ICONFIG=46  <-- CHANNEL_ID=46
    {  true, false, false, false }, // ICONFIG=47  <-- CHANNEL_ID=47
    { false, false, false,  true }, // ICONFIG=48  <-- CHANNEL_ID=48
  };

}
#endif /* clang-format on */

#endif // COLORAMPS_H
