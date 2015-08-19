#include "crystal.h"
#include <cstdlib>
using namespace std;

crystal_t::crystal_t() {

  // allow sub-objects access to crystal desc.
  symmetry.bind(this); specialk.bind(this); 
  green0.bind(this);
}

crystal_t::crystal_t(std::map<std::string,std::string>& kvp) {

  // allow sub-objects access to crystal desc.
  symmetry.bind(this); specialk.bind(this); 
  green0.bind(this);

  // configure modules
  configure(kvp);
  green0.configure(kvp);
  symmetry.configure(kvp);
  specialk.configure(kvp);
}

void crystal_t::configure(std::map<std::string,std::string>& kvp) {

  // set lattice vectors
  extract_lattvec(kvp, avec);
  set_bravais(this->avec);

  // set basis
  extract_basis(kvp, basis);
  nsites = basis.size();

  siteid.resize(nsites);
  for(int i = 0; i < nsites; i++) 
    siteid[i] = i;
 
  set_basis(this->basis,this->siteid);
}

void crystal_t::set_bravais(const mat3& _avec) {

  avec = _avec;

  // compute unit cell volume
  dble vol = dot( avec[0], cross(avec[1],avec[2]) );

  // compute reciprocal lattice vectors
  bvec[0] = cross(avec[1],avec[2])/vol;
  bvec[1] = cross(avec[2],avec[0])/vol;
  bvec[2] = cross(avec[0],avec[1])/vol;
}

void crystal_t::set_basis(const vector<vec3>& _basis, 
  const vector<int>& _siteid) {
 
  // guard against bad input
  if( _basis.size() != _siteid.size() ) { 
    printf("crystal.set_basis: basis.size() != siteid.size()\n"); 
    exit(1); 
  }

  basis = _basis;
  siteid = _siteid;

  // rewrite basis in lattice coordinates in 1st cell
  mat3 ainv = transpose(bvec);
  basisl.resize(basis.size());
  for(int i = 0; i < basis.size(); i++) {
    basisl[i] = ainv * basis[i];
    modulo_coord(basisl[i]);
  }
}

void crystal_t::find_symmetry() {
  symmetry.find();
}

void crystal_t::find_specialk(int n1, int n2, int n3) {
  specialk.generate_monkhorst_pack(n1,n2,n3);
}
