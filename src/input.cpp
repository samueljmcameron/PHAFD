#include <iostream>

#include "input.hpp"


std::vector<std::string> PHAFD::input::split_line(std::string& line)
/*
  Given an input string, split it into words (separated by whitespace) and
  return the vector of these words - removing any training comments (starting with
  '#').

 */
{
  // remove any trailing comments
  line = line.substr(0,line.find("#",0));
  
  std::vector<std::string> out;
  
  std::string tmp;
  
  
  std::size_t index;
  
  index = 0;
  while(index < line.size()) {
    
    for (; index != line.size()
	   && !isspace(line[index]); ++index)
      tmp += line[index];
    
    if (tmp != "") {
      out.push_back(tmp);
    }
    tmp = "";
    
    index += 1;
  }
  
  return out;
}


void PHAFD::input::replacePercentages(std::string &raw,int id)
/*
  Given a string which has "%" in its name, replace this with the
  the integer id.   E.g. if id = 4 and raw is

  "fname_p%.txt"

  the output would be

  "fname_p4.txt"
  
 */
{

  std::string::size_type vstart;

  while (true) {
    
    vstart = raw.find("%");
    
    if (vstart != std::string::npos) 
      raw.replace(vstart,1,std::to_string(id));
    else break;
  }
  
  return;
}
  
void PHAFD::input::convertVariables(std::string &raw,
				    std::map<std::string, std::string> const& varMap)
/*

  Given a string, replace any set of characters with the form ${var} to the
  value of var, where var and its value must be in the map varMap. E.g. for
  a var map {"name" : "Jim"}, the string

  "Hello my name is ${name}othy - well actually it's just ${name}."
  
  would convert to
  
  "Hello my name is Jimothy - well actually it's just Jim."

  
 */
{
  
  
  
  std::string::size_type vstart,vend;
  
  vend = 0;
  while (vend != raw.size()) {
    
    vstart = raw.find("${");
    
    if (vstart != std::string::npos) {
      
      vend = raw.find("}",vstart+2);
      
      if (vend != std::string::npos) {
	std::string tmp = raw.substr(vstart+2,vend-vstart-2);
	if (tmp == "")
	  throw std::runtime_error("No content ('${}') in input file.");
	bool found_key = false;
	
	for (const auto &xm: varMap) {
	  if (xm.first == tmp) {
	    found_key = true;
	    raw.erase(vstart,vend - vstart+1);
	    raw.insert(vstart,xm.second);
	  }
	}
	if (!found_key) {
	  std::string errorMessage
	    = std::string("No such command line variable ") + tmp;
	  throw std::runtime_error(errorMessage);
	}
	
	
	
      } else {
	throw::std::runtime_error("Missing '}' in input script.");
      }
      
    } else {
      vend = raw.size();
    }
  }
  
  return;
}
