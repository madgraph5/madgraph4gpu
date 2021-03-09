// See http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
#include <stdio.h>

int counters_initialise_()
{
  printf( "__CPP Initialise counters\n" );
  return(1);
}

int counters_finalise_()
{
  printf( "__CPP Finalise counters\n" );
  return(1);
}
