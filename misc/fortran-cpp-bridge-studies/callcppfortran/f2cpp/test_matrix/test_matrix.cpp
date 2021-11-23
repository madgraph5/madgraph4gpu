#include <iostream>

#define r 5
#define c 6

int main() {

  int m[r][c];

  for (int x = 1; x <= r; ++x) {
    for (int y = 1; y <= c; ++y) {
      m[x-1][y-1] = x*10 + y;
    }
  }

  for (int x = 0; x < r; ++x) {
    for (int y = 0; y < c; ++y) {
      std::cout << m[x][y] << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "--" << std::endl;

  int * is = (int*)m;
  for (int x = 0; x < r*c; ++x) {
    std::cout << is[x] << " ";
  }
  std::cout << std::endl;

  return 0;

}
