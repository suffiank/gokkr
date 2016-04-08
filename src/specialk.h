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

  /* marked for deletion
  struct star_t {
    const crystal_t* crystal;
    star_t(const crystal_t* cr) { bind(cr); }
    void bind(const crystal_t* cr) { crystal = cr; }

    vec3 k; 
    int weight() const { return symm_index.size(); }
    std::vector<int> symm_index;
    const vec3& operator[](int i) const; 
  };

  std::vector<star_t> kstar; */

  // list of k-vectors grouped in stars
  typedef std::vector<vec3> star_t;
  std::vector< star_t > kstar_kvec;
  std::vector< std::vector<int> >  kstar_symm;

public:

  void configure(std::map<std::string,std::string>& kvp);
  void bind(const crystal_t* cr); 
  void generate_monkhorst_pack(int n1, int n2, int n3);

  inline const std::vector<vec3>& operator[](int i) { return kstar_kvec[i]; }
  int numk() { return kstar_kvec.size(); };
};

#endif
