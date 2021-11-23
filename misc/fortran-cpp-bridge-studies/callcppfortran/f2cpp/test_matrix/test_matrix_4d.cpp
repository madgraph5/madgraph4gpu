#include <iostream>

using namespace std;

int gen(int part, int elem) {
  return 800 + 10*part + elem;
}

int main() {

  // AOSOA[ipagM][ipar][ip4][ieppM]
  n_evts = 2;  // ieepM
  n_part = 3;  // ipar
  n_elem = 4;  // ip4
  stride = 2;  // ipagM

  int m[n_part][n_elem];

  for (int part = 0; part < n_part; ++part) {
    for (int elem = 0; elem < n_elem; ++elem) {
      m[part][elem] = gen(part, elem);
    }
  }

  for (int part = 0; part < n_part; ++part) {
    for (int elem = 0; elem < n_elem; ++elem) {
      cout << m[part][elem] << " ";
    }
    cout << endl;
  }
  cout << endl;

  int *mp = (int*)m;
  for (int i = 0; i < part*elem; ++i) {
    cout << mp[i] << " ";
  }
  cout << endl;

}
