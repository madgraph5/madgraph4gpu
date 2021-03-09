// See http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
#include <stdio.h>

static int counters_counter = 0;

void counters_initialise_()
{
  FILE *f;
  f = fopen( "counters_log.txt", "w" );
  fprintf( f, "__CPP Initialise counters\n" );
  fclose( f );
  return;
}

void counters_start_()
{
  counters_counter++;
  return;
}

void counters_finalise_()
{
  FILE *f;
  f = fopen( "counters_log.txt", "a" );
  fprintf( f, "__CPP Finalise counters\n" );
  fprintf( f, "__CPP Total counter = %d\n", counters_counter );
  fclose( f );
  return;
}
