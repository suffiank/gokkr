#ifndef _EXTRACTKVP_H
#define _EXTRACTKVP_H

#include <string>
#include <map>
#include <vector>
#include <limits>
#include "mat3.h"

const int require_key_dbl = std::numeric_limits<double>::max();
const int require_key_int = std::numeric_limits<int>::max();

void extractkvp(std::map<std::string,std::string>& kvp, 
  const std::string& key, double &value, double defaultval); 

void extractkvp(std::map<std::string,std::string>& kvp, 
  const std::string& key, int &value, int defaultval);

void extract_lattvec(std::map<std::string,std::string>& kvp,
  mat3& avec);

void extract_basis(std::map<std::string,std::string>& kvp,
  std::vector<vec3>& basis);

#endif
