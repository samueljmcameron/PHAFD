
#ifndef PHAFD_FIXGRID_VDET_HPP
#define PHAFD_FIXGRID_VDET_HPP


#include "fix.hpp"

#include <memory>

namespace fftwArr {
  template <typename T>
  class array3D;
}

namespace PHAFD_NS {

class FixGridVdet : public Fix {
public:
  FixGridVdet(PHAFD *);

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;


  virtual void pre_final_integrate() override;
  virtual void post_final_integrate() override;
  virtual void start_of_step() override {};  
  virtual void initial_integrate() override {};
  virtual void post_force() override {};

  virtual void final_integrate() override {};

  virtual void reset_dt() override {};

  virtual void end_of_step() override {};
  
private:

  void compute_vtherm();
  void set_vdet(int, int, int) ;
  
  bool didnotintegrate;

  fftwArr::array3D<std::complex<double>> *ft_vdet_x,*ft_vdet_y, *ft_vdet_z;
  fftwArr::array3D<double> *vdet_x,*vdet_y, *vdet_z;
  
  std::vector<double> qys,qzs;

  double viscosity,temp;

};

}

#endif
