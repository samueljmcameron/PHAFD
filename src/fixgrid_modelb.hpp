
#ifndef PHAFD_FIXGRID_MODELB_HPP
#define PHAFD_FIXGRID_MODELB_HPP


#include "fix.hpp"

#include <memory>

namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {
class FixGridGradPhi;
class ConjugateNoise;


class FixGridModelB : public Fix {
public:
  FixGridModelB(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void start_of_step() override;
  
  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void pre_final_integrate() override;
  virtual void final_integrate() override;
  virtual void post_final_integrate() override;
  
  virtual void reset_dt() override;

  virtual void end_of_step() override {};
private:

  std::unique_ptr<ConjugateNoise> conjugate;

  bool didnotintegrate;
  double normalization;
  double mobility,temp,volFH,gamma;

  std::vector<double> qys,qzs;

  std::unique_ptr<FixGridGradPhi> fixgridgradphi;

  void point_update(int,int,int);

  void origin_update();

  fftwArr::array3D<std::complex<double>> *ft_phi;
  fftwArr::array3D<std::complex<double>> *ft_chempot;

  
  
};

}

#endif
