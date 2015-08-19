#ifndef _SPECIALK_H
#define _SPECIALK_H

#include <vector>
#include "symmetry.h"

// note: specialk[i][j] returns jth k-vector of ith star
// note: this requires a matrix operation for each access
// note: also requires binding crystal* to each k-point (redundant)

class specialk_t {

  const crystal_t* crystal;

  void fold_specialk(std::vector<vec3>& specialk);
  void reduce_to_bzone(vec3& kpoint);

  struct star_t {
    const crystal_t* crystal;
    star_t(const crystal_t* cr) { bind(cr); }
    void bind(const crystal_t* cr) { crystal = cr; }

    vec3 k; std::vector<int> symm_index;
    vec3 operator[](int i); 
  };

  std::vector<star_t> kstar;

public:

  void configure(std::map<std::string,std::string>& kvp);
  void bind(const crystal_t* cr); 
  void generate_monkhorst_pack(int n1, int n2, int n3);

  inline const star_t& operator[](int i) { return kstar[i]; }
};

#endif
