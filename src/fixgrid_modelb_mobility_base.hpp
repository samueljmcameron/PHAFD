
#ifndef PHAFD_FIXGRID_MODELB_MOBILITY_BASE_HPP
#define PHAFD_FIXGRID_MODELB_MOBILITY_BASE_HPP


#include "fix.hpp"
#include <fftw3-mpi.h>

#include <memory>
#include <complex>

namespace fftwArr {
  template<typename>
  class array3D;
}

namespace PHAFD_NS {
class FixGridGradient;
class ConjugateNoiseNoQ;
class FixGridDivergence;


class FixGridModelBMobilityBase : public Fix {
public:
  FixGridModelBMobilityBase(PHAFD *);
  ~FixGridModelBMobilityBase();

  virtual void init(const std::vector<std::string> &) override;
  
  virtual void setup() override;

  virtual void start_of_step() override;
  
  virtual void initial_integrate() override {};
  virtual void post_force() override {};
  virtual void pre_final_integrate() override;
  virtual void final_integrate() override;
  virtual void post_final_integrate(bool invert_fft = true) override;
  
  virtual void reset_dt() override;

  virtual void end_of_step() override {};
private:

  double temp; // temperature (necessary for stochastic drift part)
  std::array<std::unique_ptr<ConjugateNoiseNoQ>,3> conjugate_noise;

  std::unique_ptr<FixGridGradient> gradfix;
  std::unique_ptr<FixGridDivergence> divfix;

  
  double normalization;
  double mobility;

  void point_update(int,int,int);

  void origin_update();
  std::array<std::unique_ptr<fftwArr::array3D<double>>,3> flux;
  std::array<std::unique_ptr<fftwArr::array3D<std::complex<double>>>,
	     3> ft_flux;

  std::array<fftw_plan,3> forward_flux;
  void compute_stochastic_drift();
  void compute_usual_drift();
  void add_noise_to_phi();
  bool plan_set;
  std::array<fftwArr::array3D<std::complex<double>> *,3> ft_rnoises;
  std::array<std::unique_ptr<fftwArr::array3D<double>> ,3> rnoises;

  //std::unique_ptr<fftwArr::array3D<double>> noise;
  // std::unique_ptr<fftwArr::array3D<std::complex<double>>> ft_noise;

  //fftw_plan forward_noise;
  
  std::array<fftw_plan,3> backward_rnoises;
  double inv_vol_element;
  fftw_plan forward_mobility_deriv, forward_sqrt_mobility;
  std::unique_ptr<fftwArr::array3D<std::complex<double>>> ft_sqrt_mobility;
  std::unique_ptr<fftwArr::array3D<std::complex<double>>> ft_mobility_deriv;

  
protected:
  std::unique_ptr<fftwArr::array3D<double>> sqrt_mobility;

  std::unique_ptr<fftwArr::array3D<double>> mobility_deriv;

  std::array<std::unique_ptr<fftwArr::array3D<double>>,3> grad_sqrt_mobility;

  virtual void calculate_sqrt_mobility() = 0;
  virtual void calculate_mobility_deriv() = 0;



  
};

}

#endif
