#ifndef PHAFD_DOMAIN_HPP
#define PHAFD_DOMAIN_HPP


#include "atom.hpp"
#include <array>
#include <Eigen/Core>


#include "ps_pde/domain.hpp"
#include "ps_pde/grid.hpp"

namespace PHAFD {
class Domain : public psPDE::Domain
{
public:
  Domain(int,int, std::string);
  
  void partition(const psPDE::Grid &);
  void pbc (Atom &) const;
  int set_image() const;


  void map(Eigen::Ref<Eigen::Vector3d>,
  	   const Eigen::Ref<const Eigen::Vector3d> &,int ) const;

private:
  void unmap(Eigen::Ref<Eigen::Vector3d>,
  	     const Eigen::Ref<const Eigen::Vector3d> &,int image) const;

};
}
#endif
