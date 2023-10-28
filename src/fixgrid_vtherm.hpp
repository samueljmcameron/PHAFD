
#ifndef PHAFD_FIXGRID_VTHERM_HPP
#define PHAFD_FIXGRID_VTHERM_HPP


#include "fix.hpp"

#include <memory>

namespace fftwArr {
  template <typename T>
  class array3D;
}

namespace PHAFD_NS {

class ConjugateNoise;

class FixGridVtherm : public Fix {
public:
  FixGridVtherm(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void start_of_step() override;
  
  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  
  virtual void reset_dt() override;

  virtual void end_of_step() override {};
  
private:

  void compute_vtherm();
  void set_vtherm(int, int, int) ;
  
  std::array<std::unique_ptr<ConjugateNoise>,3> conjugate_vnoise;

  bool didnotintegrate;

  fftwArr::array3D<std::complex<double>> *ft_Znoise_x,*ft_Znoise_y, *ft_Znoise_z;
  fftwArr::array3D<std::complex<double>> *ft_vtherm_x,*ft_vtherm_y, *ft_vtherm_z;
  
  std::vector<double> qys,qzs;

  double viscosity,temp;

  
};

}

#endif
