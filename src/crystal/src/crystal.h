#ifndef _CRYSTAL_H
#define _CRYSTAL_H

#include <vector>
#include <map>
#include <string>
#include "../../util/src/mat3.h"
#include "../../util/src/mathfn.h"
#include "../../util/src/extractkvp.h"

#include "symmetry.h"
#include "specialk.h"
#include "../../green0/src/green0.h"
#include "../../atom/src/atom.h"

struct crystal_t {

  mat3 avec, bvec;
  int nsites;
  std::vector<vec3> basis;
  std::vector<vec3> basisl;
  std::vector<int> siteid;
  std::vector<atom_t> atom;
  dble vol, bzvol;

  crystal_t();
  crystal_t(std::map<std::string,std::string>& kvp);
  void configure(std::map<std::string,std::string>& kvp);

  void set_bravais(const mat3& _avec);
  void set_basis(const std::vector<vec3>& _basis,
    const std::vector<int>& _siteid);

  symmetry_t symmetry;
  void find_symmetry();

  specialk_t specialk;
  void find_specialk(int n1, int n2, int n3);

  green0_t green0;

  // main methods
  cplx lloydN(cplx E);
};

#endif
