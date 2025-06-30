#ifndef PHAFD_CONJUGATE_NOISE_NO_Q_HPP
#define PHAFD_CONJUGATE_NOISE_NO_Q_HPP

#include <random>

#include "conjugate_noise.hpp"

#include <complex>
#include <vector>
#include <string>


namespace PHAFD_NS {


class ConjugateNoiseNoQ : public ConjugateNoise {
public:
  ConjugateNoiseNoQ(PHAFD *);

private:
  virtual void complex_update(int,int,int) override;
  virtual void real_update(int,int,int) override;

};

}

#endif
