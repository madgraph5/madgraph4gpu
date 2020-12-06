#include "alsimple.h"

namespace alsimple {

std::random_device rd;

randGenerator::randGenerator() : gen(rd()), dis_double(0.0,1.0), dis_float(0.0f,1.0f) { }

randStatus_t createGeneratorHost( randGenerator_t *gen, alsrandType_t type ) {
  *gen = new randGenerator;
  return RANDOK;
}

randStatus_t createGenerator( randGenerator_t *gen, alsrandType_t type ) {
  // not supported
  return RANDBAD;
}

randStatus_t randSetPseudoRandomGeneratorSeed( randGenerator_t gen, unsigned long long seed ) {
  std::seed_seq seq{seed, seed>>32};
  gen->gen.seed(seq);
  return RANDOK;
}

randStatus_t randDestroyGenerator( randGenerator_t gen ) {
  delete gen;
  return RANDOK;
}

randStatus_t randGenerateUniformDouble( randGenerator_t gen, double *out, size_t nsize) {
  for(size_t i=0;i<nsize;++i) {
    out[i] = gen->dis_double(gen->gen);
  }
  return RANDOK;
}

randStatus_t randGenerateUniform( randGenerator_t gen, float *out, size_t nsize) {
  for(size_t i=0;i<nsize;++i) {
    out[i] = gen->dis_float(gen->gen);
  }
  return RANDOK;
}

}
