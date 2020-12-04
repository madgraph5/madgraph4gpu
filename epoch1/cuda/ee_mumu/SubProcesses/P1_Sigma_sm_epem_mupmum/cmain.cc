int check( int argc, char **argv ); // from check.cc (compiled with g++)
int gcheck( int argc, char **argv ); // from gcheck.cu (compiled with nvcc)

// This is built with g++ and linked with objects compiled with g++
int main( int argc, char **argv )
{
  return check( argc, argv );
}
