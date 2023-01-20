#ifndef PHAFD_INPUT_HPP
#define PHAFD_INPUT_HPP

#include <string>
#include <vector>
#include <map>


namespace PHAFD {
  namespace input {
  std::vector<std::string> split_line(std::string&);
  void replacePercentages(std::string &, int);

  void convertVariables(std::string &,
			std::map<std::string, std::string> const&);
  }
  
};
#endif
