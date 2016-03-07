#ifndef _SYMMETRY_H
#define _SYMMETRY_H

#include <string>
#include <map>
#include <vector>
#include "../../util/src/mat3.h"
#include "../../util/src/mathfn.h"
class crystal_t;

class symmetry_t {

  static bool init_std_symm;
  std::vector<mat3> cube_symmetry;
  std::vector<mat3> hex_symmetry;

  const crystal_t* crystal;

  std::vector<mat3> crystal_symmetry;
  int nsymm;

  std::string system, basegroup;
  std::vector<int> subgroup;
  int suborder;

  void fill_cube_isometries(std::vector<mat3>& rotmat);
  void fill_hex_isometries(std::vector<mat3>& rotmat);
  void fill_rotation_axis(mat3& rotmat, dble theta, const vec3& axis);
  void fill_period(std::vector<mat3>& rotmat, int& n);

  void find_bravais_symmetry();
  void find_crystal_symmetry();

public:

  symmetry_t();

  void configure(std::map<std::string,std::string>& kvp);
  void bind(const crystal_t* cryst) { crystal = cryst; }
  void find() { find_crystal_symmetry(); }
  void make_sakuraiD();

  inline mat3 operator[](int i) const { return crystal_symmetry[i]; }
  inline int size() const { return crystal_symmetry.size(); }
  
  std::vector< std::vector<cplx> > sakuraiD;
};

#endif
