#ifndef PHAFD_DOMAIN_HPP

#define PHAFD_DOMAIN_HPP


#include <array>
#include <Eigen/Core>

#include "pointers.hpp"

namespace PHAFD_NS {

  
class Domain : protected Pointers
{
public:
  Domain(PHAFD *);
  ~Domain();
  
  void set_box(const std::vector<std::string> &);
  void set_subbox();
  void pbc () const;
  int set_image() const;

  bool boxset, subboxset;
  
  // addresses of data from psPDE arrays (of size 3)
  std::array<double,3> period,boxlo,boxhi;
  std::array<double,3> sublo,subhi;

  void map(Eigen::Ref<Eigen::Vector3d>,
  	   const Eigen::Ref<const Eigen::Vector3d> &,int ) const;

  double dqx();
  double dqy();
  double dqz();
  
private:

  void unmap(Eigen::Ref<Eigen::Vector3d>,
  	     const Eigen::Ref<const Eigen::Vector3d> &,int image) const;

};
}
#endif
