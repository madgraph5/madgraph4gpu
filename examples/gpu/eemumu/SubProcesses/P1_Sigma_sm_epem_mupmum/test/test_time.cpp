#include <iostream>
#include "timer.h"
#include <unistd.h>
#include <iomanip>


int main() {
 
  Timer<std::chrono::high_resolution_clock> t;

  t.Info();

  t.Start();
  usleep(1000);
  float f = t.GetDuration();

  std::cout << f << std::endl;
 

  return 0;
}
