#ifndef PHAFD_CONJUGATE_VOLFRAC_HPP
#define PHAFD_CONJUGATE_VOLFRAC_HPP

#include <random>

#include "conjugate.hpp"

#include <complex>
#include <vector>
#include <string>

namespace fftwArr {
  template<typename>
  class array3D;

}

namespace PHAFD_NS {


class ConjugateVolFrac : public Conjugate {
public:
  ConjugateVolFrac(PHAFD *);

  virtual void readCoeffs(const std::vector<std::string> &) override;
  

  virtual void reset_dt(double) override;
private:


  std::mt19937 gen;

  std::uniform_real_distribution<double> real_dist;

  double normalization;

  double mobility,temp,volFH,gamma,dt;

  int seed;

  bool seed_flag;

  std::complex<double> noise;

  fftwArr::array3D<std::complex<double>> *ft_phi;
  fftwArr::array3D<std::complex<double>> *ft_chempot;


  
  
  virtual void complex_update(int,int,int) override;
  virtual void real_update(int,int,int) override;
  virtual void origin_update() override;

};

}

#endif
