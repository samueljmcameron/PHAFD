#ifndef PHAFD_CONJUGATE_NOISE_HPP
#define PHAFD_CONJUGATE_NOISE_HPP

#include <random>

#include "conjugate.hpp"

#include <complex>
#include <vector>
#include <string>


namespace PHAFD_NS {


class ConjugateNoise : public Conjugate {
public:
  ConjugateNoise(PHAFD *);

  virtual void readCoeffs(const std::vector<std::string> &) override;
  

  virtual void reset_dt(double) override;

  double damping,temp;
private:


  std::mt19937 gen;

  std::uniform_real_distribution<double> real_dist;

  double dt;

  int seed;
  bool seed_flag;
  
  virtual void complex_update(int,int,int) override;
  virtual void real_update(int,int,int) override;
  virtual void origin_update() override;

};

}

#endif
