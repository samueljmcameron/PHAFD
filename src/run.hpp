#ifndef RUN_HPP
#define RUN_HPP

#include <vector>
#include <string>
#include "globalparams.hpp"
#include "ps_pde/solutionparams.hpp"

void run(GlobalParams,psPDE::SolutionParams,const std::vector<std::string> &,
	 const std::vector<std::vector<std::string>> &,
	 std::vector<std::vector<double>> &);




#endif
