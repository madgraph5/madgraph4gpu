// See http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
#include <stdio.h>

static int counters_counter1 = 0;
static int counters_counter2 = 0;

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
  counters_counter1++;
  return;
}

void counters_end_()
{
  counters_counter2++;
  return;
}

void counters_finalise_()
{
  FILE *f;
  f = fopen( "counters_log.txt", "a" );
  fprintf( f, "__CPP Finalise counters\n" );
  fprintf( f, "__CPP Total counter1 = %d\n", counters_counter1 );
  fprintf( f, "__CPP Total counter2 = %d\n", counters_counter2 );
  fclose( f );
  return;
}
