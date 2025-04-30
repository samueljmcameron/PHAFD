
#ifndef PHAFD_FIXGRID_MODELH_HPP
#define PHAFD_FIXGRID_MODELH_HPP


#include "fix.hpp"

#include <memory>

namespace fftwArr {
  template <typename T>
  class array3D;
}

namespace PHAFD_NS {

class ConjugateNoise;
class FixGridGradPhi;
class FixGridVtherm;

class FixGridModelH : public Fix {
public:
  FixGridModelH(PHAFD *);

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

  std::unique_ptr<ConjugateNoise> conjugate_phinoise;
  std::array<std::unique_ptr<ConjugateNoise>,3> conjugate_vnoise;
  std::unique_ptr<FixGridGradPhi> fixgridgradphi;
  std::unique_ptr<FixGridVtherm> fixgridvtherm;
  bool didnotintegrate;

  double *qys,qzs;


  fftwArr::array3D<std::complex<double>> *ft_Znoise_x,*ft_Znoise_y, *ft_Znoise_z;

  
};

}

#endif
