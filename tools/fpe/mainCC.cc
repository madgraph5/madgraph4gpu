extern "C" void testfpes_( const bool* pClearFPEs );

int main()
{
  const bool clearFPEs = true;
  //const bool clearFPEs = false;
  testfpes_( &clearFPEs );
  return 0;
}
