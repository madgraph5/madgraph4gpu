int check( int argc, char **argv ); // from check.cc (compiled with g++)
int gcheck( int argc, char **argv ); // from gcheck.cu (compiled with nvcc)

// This is compiled with g++ and linked with objects compiled with nvcc or g++
int main( int argc, char **argv )
{
  if ( gcheck( argc, argv ) != 0 ) return 1;
  if ( check( argc, argv ) != 0 ) return 2;
  return 0;
}
