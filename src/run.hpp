#ifndef RUN_HPP
#define RUN_HPP

#include <vector>
#include <string>
#include "ps_pde/globalparams.hpp"
#include "ps_pde/solutionparams.hpp"

void run(psPDE::GlobalParams,psPDE::SolutionParams,const std::vector<std::string> &,
	 const std::vector<std::vector<std::string>> &,
	 std::vector<std::vector<double>> &,std::vector<double> &,
	 std::vector<double> &);




#endif
