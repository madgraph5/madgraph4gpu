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
  __device__ constexpr int channel2iconfig[96] = { // note: a trailing comma in the initializer list is allowed
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
     49, // CHANNEL_ID=49  i.e. DIAGRAM=49  --> ICONFIG=49
     50, // CHANNEL_ID=50  i.e. DIAGRAM=50  --> ICONFIG=50
     51, // CHANNEL_ID=51  i.e. DIAGRAM=51  --> ICONFIG=51
     52, // CHANNEL_ID=52  i.e. DIAGRAM=52  --> ICONFIG=52
     53, // CHANNEL_ID=53  i.e. DIAGRAM=53  --> ICONFIG=53
     54, // CHANNEL_ID=54  i.e. DIAGRAM=54  --> ICONFIG=54
     55, // CHANNEL_ID=55  i.e. DIAGRAM=55  --> ICONFIG=55
     56, // CHANNEL_ID=56  i.e. DIAGRAM=56  --> ICONFIG=56
     57, // CHANNEL_ID=57  i.e. DIAGRAM=57  --> ICONFIG=57
     58, // CHANNEL_ID=58  i.e. DIAGRAM=58  --> ICONFIG=58
     59, // CHANNEL_ID=59  i.e. DIAGRAM=59  --> ICONFIG=59
     60, // CHANNEL_ID=60  i.e. DIAGRAM=60  --> ICONFIG=60
     61, // CHANNEL_ID=61  i.e. DIAGRAM=61  --> ICONFIG=61
     62, // CHANNEL_ID=62  i.e. DIAGRAM=62  --> ICONFIG=62
     63, // CHANNEL_ID=63  i.e. DIAGRAM=63  --> ICONFIG=63
     64, // CHANNEL_ID=64  i.e. DIAGRAM=64  --> ICONFIG=64
     65, // CHANNEL_ID=65  i.e. DIAGRAM=65  --> ICONFIG=65
     66, // CHANNEL_ID=66  i.e. DIAGRAM=66  --> ICONFIG=66
     67, // CHANNEL_ID=67  i.e. DIAGRAM=67  --> ICONFIG=67
     68, // CHANNEL_ID=68  i.e. DIAGRAM=68  --> ICONFIG=68
     69, // CHANNEL_ID=69  i.e. DIAGRAM=69  --> ICONFIG=69
     70, // CHANNEL_ID=70  i.e. DIAGRAM=70  --> ICONFIG=70
     71, // CHANNEL_ID=71  i.e. DIAGRAM=71  --> ICONFIG=71
     72, // CHANNEL_ID=72  i.e. DIAGRAM=72  --> ICONFIG=72
     73, // CHANNEL_ID=73  i.e. DIAGRAM=73  --> ICONFIG=73
     74, // CHANNEL_ID=74  i.e. DIAGRAM=74  --> ICONFIG=74
     75, // CHANNEL_ID=75  i.e. DIAGRAM=75  --> ICONFIG=75
     76, // CHANNEL_ID=76  i.e. DIAGRAM=76  --> ICONFIG=76
     77, // CHANNEL_ID=77  i.e. DIAGRAM=77  --> ICONFIG=77
     78, // CHANNEL_ID=78  i.e. DIAGRAM=78  --> ICONFIG=78
     79, // CHANNEL_ID=79  i.e. DIAGRAM=79  --> ICONFIG=79
     80, // CHANNEL_ID=80  i.e. DIAGRAM=80  --> ICONFIG=80
     81, // CHANNEL_ID=81  i.e. DIAGRAM=81  --> ICONFIG=81
     82, // CHANNEL_ID=82  i.e. DIAGRAM=82  --> ICONFIG=82
     83, // CHANNEL_ID=83  i.e. DIAGRAM=83  --> ICONFIG=83
     84, // CHANNEL_ID=84  i.e. DIAGRAM=84  --> ICONFIG=84
     85, // CHANNEL_ID=85  i.e. DIAGRAM=85  --> ICONFIG=85
     86, // CHANNEL_ID=86  i.e. DIAGRAM=86  --> ICONFIG=86
     87, // CHANNEL_ID=87  i.e. DIAGRAM=87  --> ICONFIG=87
     88, // CHANNEL_ID=88  i.e. DIAGRAM=88  --> ICONFIG=88
     89, // CHANNEL_ID=89  i.e. DIAGRAM=89  --> ICONFIG=89
     90, // CHANNEL_ID=90  i.e. DIAGRAM=90  --> ICONFIG=90
     91, // CHANNEL_ID=91  i.e. DIAGRAM=91  --> ICONFIG=91
     92, // CHANNEL_ID=92  i.e. DIAGRAM=92  --> ICONFIG=92
     93, // CHANNEL_ID=93  i.e. DIAGRAM=93  --> ICONFIG=93
     94, // CHANNEL_ID=94  i.e. DIAGRAM=94  --> ICONFIG=94
     95, // CHANNEL_ID=95  i.e. DIAGRAM=95  --> ICONFIG=95
     96, // CHANNEL_ID=96  i.e. DIAGRAM=96  --> ICONFIG=96
  };

  // Map iconfig to the mask of allowed colors (e.g. "colormask = icolamp[iconfig - 1]": input index uses C indexing)
  // This array has N_config <= N_diagrams elements
  __device__ constexpr bool icolamp[96][6] = { // note: a trailing comma in the initializer list is allowed
    {  true, false, false, false, false, false }, // ICONFIG=1   <-- CHANNEL_ID=1
    { false,  true, false, false, false, false }, // ICONFIG=2   <-- CHANNEL_ID=2
    {  true, false, false, false, false, false }, // ICONFIG=3   <-- CHANNEL_ID=3
    { false,  true, false, false, false, false }, // ICONFIG=4   <-- CHANNEL_ID=4
    {  true,  true, false, false, false, false }, // ICONFIG=5   <-- CHANNEL_ID=5
    {  true,  true, false, false, false, false }, // ICONFIG=6   <-- CHANNEL_ID=6
    {  true,  true, false, false, false, false }, // ICONFIG=7   <-- CHANNEL_ID=7
    {  true,  true, false, false, false, false }, // ICONFIG=8   <-- CHANNEL_ID=8
    { false,  true, false, false, false, false }, // ICONFIG=9   <-- CHANNEL_ID=9
    { false,  true, false, false, false, false }, // ICONFIG=10  <-- CHANNEL_ID=10
    { false,  true, false, false, false, false }, // ICONFIG=11  <-- CHANNEL_ID=11
    { false,  true, false, false, false, false }, // ICONFIG=12  <-- CHANNEL_ID=12
    {  true, false, false, false, false, false }, // ICONFIG=13  <-- CHANNEL_ID=13
    {  true, false, false, false, false, false }, // ICONFIG=14  <-- CHANNEL_ID=14
    {  true, false, false, false, false, false }, // ICONFIG=15  <-- CHANNEL_ID=15
    {  true, false, false, false, false, false }, // ICONFIG=16  <-- CHANNEL_ID=16
    { false, false, false, false,  true,  true }, // ICONFIG=17  <-- CHANNEL_ID=17
    { false, false, false, false,  true,  true }, // ICONFIG=18  <-- CHANNEL_ID=18
    { false, false, false, false,  true,  true }, // ICONFIG=19  <-- CHANNEL_ID=19
    { false, false, false, false,  true,  true }, // ICONFIG=20  <-- CHANNEL_ID=20
    {  true, false,  true, false, false, false }, // ICONFIG=21  <-- CHANNEL_ID=21
    {  true, false,  true, false,  true,  true }, // ICONFIG=22  <-- CHANNEL_ID=22
    {  true, false,  true, false,  true,  true }, // ICONFIG=23  <-- CHANNEL_ID=23
    { false, false, false, false,  true,  true }, // ICONFIG=24  <-- CHANNEL_ID=24
    {  true, false,  true, false, false, false }, // ICONFIG=25  <-- CHANNEL_ID=25
    {  true, false,  true, false,  true,  true }, // ICONFIG=26  <-- CHANNEL_ID=26
    {  true, false,  true, false,  true,  true }, // ICONFIG=27  <-- CHANNEL_ID=27
    { false, false, false, false,  true,  true }, // ICONFIG=28  <-- CHANNEL_ID=28
    {  true, false,  true, false, false, false }, // ICONFIG=29  <-- CHANNEL_ID=29
    {  true, false,  true, false, false, false }, // ICONFIG=30  <-- CHANNEL_ID=30
    {  true, false,  true, false, false, false }, // ICONFIG=31  <-- CHANNEL_ID=31
    {  true, false,  true, false, false, false }, // ICONFIG=32  <-- CHANNEL_ID=32
    { false, false,  true,  true, false, false }, // ICONFIG=33  <-- CHANNEL_ID=33
    { false, false,  true,  true, false, false }, // ICONFIG=34  <-- CHANNEL_ID=34
    { false, false,  true,  true, false, false }, // ICONFIG=35  <-- CHANNEL_ID=35
    { false, false,  true,  true, false, false }, // ICONFIG=36  <-- CHANNEL_ID=36
    { false,  true, false, false,  true, false }, // ICONFIG=37  <-- CHANNEL_ID=37
    { false,  true,  true,  true,  true, false }, // ICONFIG=38  <-- CHANNEL_ID=38
    { false,  true,  true,  true,  true, false }, // ICONFIG=39  <-- CHANNEL_ID=39
    { false, false,  true,  true, false, false }, // ICONFIG=40  <-- CHANNEL_ID=40
    { false,  true, false, false,  true, false }, // ICONFIG=41  <-- CHANNEL_ID=41
    { false,  true,  true,  true,  true, false }, // ICONFIG=42  <-- CHANNEL_ID=42
    { false,  true,  true,  true,  true, false }, // ICONFIG=43  <-- CHANNEL_ID=43
    { false, false,  true,  true, false, false }, // ICONFIG=44  <-- CHANNEL_ID=44
    { false,  true, false, false,  true, false }, // ICONFIG=45  <-- CHANNEL_ID=45
    { false,  true, false, false,  true, false }, // ICONFIG=46  <-- CHANNEL_ID=46
    { false,  true, false, false,  true, false }, // ICONFIG=47  <-- CHANNEL_ID=47
    { false,  true, false, false,  true, false }, // ICONFIG=48  <-- CHANNEL_ID=48
    { false, false, false,  true, false, false }, // ICONFIG=49  <-- CHANNEL_ID=49
    { false, false, false,  true, false, false }, // ICONFIG=50  <-- CHANNEL_ID=50
    { false, false, false,  true, false, false }, // ICONFIG=51  <-- CHANNEL_ID=51
    { false, false, false,  true, false, false }, // ICONFIG=52  <-- CHANNEL_ID=52
    { false, false, false, false, false,  true }, // ICONFIG=53  <-- CHANNEL_ID=53
    { false, false, false, false, false,  true }, // ICONFIG=54  <-- CHANNEL_ID=54
    { false, false, false, false, false,  true }, // ICONFIG=55  <-- CHANNEL_ID=55
    { false, false, false, false, false,  true }, // ICONFIG=56  <-- CHANNEL_ID=56
    { false, false, false, false, false,  true }, // ICONFIG=57  <-- CHANNEL_ID=57
    { false, false, false,  true, false, false }, // ICONFIG=58  <-- CHANNEL_ID=58
    { false, false, false, false, false,  true }, // ICONFIG=59  <-- CHANNEL_ID=59
    { false, false, false,  true, false, false }, // ICONFIG=60  <-- CHANNEL_ID=60
    { false, false, false,  true, false,  true }, // ICONFIG=61  <-- CHANNEL_ID=61
    { false, false, false,  true, false,  true }, // ICONFIG=62  <-- CHANNEL_ID=62
    { false, false, false,  true, false,  true }, // ICONFIG=63  <-- CHANNEL_ID=63
    { false, false, false,  true, false,  true }, // ICONFIG=64  <-- CHANNEL_ID=64
    { false, false,  true, false, false, false }, // ICONFIG=65  <-- CHANNEL_ID=65
    { false, false, false,  true, false, false }, // ICONFIG=66  <-- CHANNEL_ID=66
    { false, false,  true, false, false, false }, // ICONFIG=67  <-- CHANNEL_ID=67
    { false, false, false,  true, false, false }, // ICONFIG=68  <-- CHANNEL_ID=68
    { false, false,  true, false, false, false }, // ICONFIG=69  <-- CHANNEL_ID=69
    { false, false,  true, false, false, false }, // ICONFIG=70  <-- CHANNEL_ID=70
    { false, false,  true, false, false, false }, // ICONFIG=71  <-- CHANNEL_ID=71
    { false, false,  true, false, false, false }, // ICONFIG=72  <-- CHANNEL_ID=72
    { false, false, false, false,  true, false }, // ICONFIG=73  <-- CHANNEL_ID=73
    { false, false, false, false, false,  true }, // ICONFIG=74  <-- CHANNEL_ID=74
    { false, false, false, false,  true, false }, // ICONFIG=75  <-- CHANNEL_ID=75
    { false, false, false, false, false,  true }, // ICONFIG=76  <-- CHANNEL_ID=76
    { false, false, false, false,  true, false }, // ICONFIG=77  <-- CHANNEL_ID=77
    { false, false, false, false,  true, false }, // ICONFIG=78  <-- CHANNEL_ID=78
    { false, false, false, false,  true, false }, // ICONFIG=79  <-- CHANNEL_ID=79
    { false, false, false, false,  true, false }, // ICONFIG=80  <-- CHANNEL_ID=80
    {  true,  true, false,  true, false,  true }, // ICONFIG=81  <-- CHANNEL_ID=81
    {  true,  true, false,  true, false,  true }, // ICONFIG=82  <-- CHANNEL_ID=82
    {  true,  true, false, false, false, false }, // ICONFIG=83  <-- CHANNEL_ID=83
    { false, false, false,  true, false,  true }, // ICONFIG=84  <-- CHANNEL_ID=84
    {  true,  true, false,  true, false,  true }, // ICONFIG=85  <-- CHANNEL_ID=85
    {  true,  true, false,  true, false,  true }, // ICONFIG=86  <-- CHANNEL_ID=86
    {  true,  true, false, false, false, false }, // ICONFIG=87  <-- CHANNEL_ID=87
    { false, false, false,  true, false,  true }, // ICONFIG=88  <-- CHANNEL_ID=88
    { false, false, false, false,  true, false }, // ICONFIG=89  <-- CHANNEL_ID=89
    { false,  true, false, false, false, false }, // ICONFIG=90  <-- CHANNEL_ID=90
    { false, false, false, false,  true, false }, // ICONFIG=91  <-- CHANNEL_ID=91
    { false,  true, false, false, false, false }, // ICONFIG=92  <-- CHANNEL_ID=92
    { false, false,  true, false, false, false }, // ICONFIG=93  <-- CHANNEL_ID=93
    {  true, false, false, false, false, false }, // ICONFIG=94  <-- CHANNEL_ID=94
    { false, false,  true, false, false, false }, // ICONFIG=95  <-- CHANNEL_ID=95
    {  true, false, false, false, false, false }, // ICONFIG=96  <-- CHANNEL_ID=96
  };

}
#endif /* clang-format on */

#endif // COLORAMPS_H
