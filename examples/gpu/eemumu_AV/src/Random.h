class Random {
public:
  double ranmar();
  void rmarin(int ij, int kl);

  double ranmar2() // avoid values <1e-16 as in rn()
  {
    while (true) {
      double ran = ranmar();
      if (ran > 1e-16) return ran;
    }
  }

private:
  double ranu[98];
  double ranc, rancd, rancm;
  int iranmr, jranmr;
};

double rn(int idummy);
