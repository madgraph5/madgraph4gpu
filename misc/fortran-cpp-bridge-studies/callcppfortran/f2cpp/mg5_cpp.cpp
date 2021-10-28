#include <complex>
#include <iostream>

using namespace std;

void foo_mi2(int* m, int r, int c) {
  cout << endl;

  for (int i = 0; i < r*c; ++i) {
    cout << m[i] << " ";
  }
  cout << endl << endl;

  int (*m2)[r] = (int (*)[r])m;
  for (int i = 0; i < c; ++i) {
    for (int j = 0; j < r; ++j) {
      cout << m2[i][j] << " ";
    }
    cout << endl;
  }
}

extern "C" {

  void foo_i_(int* n) {
    for (int i = 0; i < 3; ++i)
      n[i] *= -1;
  }

  void foo_f_(float* n) {
    for (int i = 0; i < 3; ++i)
      n[i] *= 2;
  }

  void foo_d_(double* n) {
    for (int i = 0; i < 3; ++i)
      n[i] /= 2;
  }

  void foo_cf_(std::complex<float>* n) {
    for (int i = 0; i < 3; ++i)
      n[i] += 2;
  }

  void foo_cd_(std::complex<double>* n) {
    for (int i = 0; i < 3; ++i)
      n[i] -= 2;
  }

  void foo_mi_(int* m, int *r, int *c) {
    foo_mi2(m,*r,*c);
  }

}
