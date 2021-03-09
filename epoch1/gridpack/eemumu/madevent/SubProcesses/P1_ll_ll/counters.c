// See http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
#include <stdio.h>

int counters_initialise_()
{
  FILE *f;
  f = fopen( "counters_log.txt", "w" );
  fprintf( f, "__CPP Initialise counters\n" );
  fclose( f );
  return( 1 );
}

int counters_finalise_()
{
  FILE *f;
  f = fopen( "counters_log.txt", "a" );
  fprintf( f, "__CPP Finalise counters\n" );
  fclose( f );
  return( 1 );
}
