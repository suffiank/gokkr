/* written by Suffian Khan Aug 2014 */

#include "green0.h"
#include "gsl/gsl_sort.h"
using namespace std;

// generate all lattice vectors within sphere of radius 'rad'
void green0_t::genlattice(vector<vec3> *latt, 
  const mat3& avec, const double& rad) {

  // need tangent planes to find bounds that
  // determine lattice circumscribing sphere

  // tangent planes to sphere given by
  //   n1 a1 + u a2 + v a3 = hat(a2 x a3) R
  // solve for n1 by multiplying by (a2 x a3)
  // solution is mag(b1)/2pi for reciprocal vector b1
  double vol = fabs(dot(avec[0],cross(avec[1],avec[2])));
  int n1 = ceil(mag(cross(avec[1],avec[2]))*rad/vol);
  int n2 = ceil(mag(cross(avec[2],avec[0]))*rad/vol);
  int n3 = ceil(mag(cross(avec[0],avec[1]))*rad/vol);

  // for every vector in the circumscribing lattice 
  vector<vec3> ulatt; vector<double> vlen;
  for(int i1 = -n1; i1 <= n1; i1++)
  for(int i2 = -n2; i2 <= n2; i2++)
  for(int i3 = -n3; i3 <= n3; i3++) {
    vec3 v = i1*avec[0] + i2*avec[1] + i3*avec[2];
    
    // keep those vectors that are inside sphere
    double magv = mag(v);
    if( mag(v) <= rad ) {
      ulatt.push_back(v);
      vlen.push_back(magv);
    }
  }

  // todo: replace with an in-place sort using STL?
  // find the permutation to sort vectors by length
  vector<size_t> perm(vlen.size());
  gsl_sort_index(perm.data(),vlen.data(),1,vlen.size());

  // store final vectors according to sorted permutation
  latt->clear(); latt->resize(vlen.size());
  for(int i = 0; i < perm.size(); i++)
    (*latt)[i] = ulatt[perm[i]];

}


