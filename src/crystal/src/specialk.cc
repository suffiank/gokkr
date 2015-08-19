#include "crystal.h"
using namespace std;

static const dble tol = 1.e-06;

void specialk_t::configure(std::map<std::string,std::string>& kvp) {

}

void specialk_t::bind(const crystal_t* cr) {
  crystal = cr;
  for(int i = 0; i < kstar.size(); i++)
    kstar[i].bind(cr);
}

inline vec3 specialk_t::star_t::operator[](int i) { 
  return crystal->symmetry[ symm_index[i] ] * k;
}
 
void specialk_t::generate_monkhorst_pack(int n1, int n2, int n3) {

    // generate specialk mesh within brillouin zone
    int n = 0; vector<vec3> specialk;
    for(int i1 = 0; i1 < n1; i1++)
    for(int i2 = 0; i2 < n2; i2++)
    for(int i3 = 0; i3 < n3; i3++) {

      dble u1 = dble(-1+n1-2*i1)/dble(2*n1);
      dble u2 = dble(-1+n2-2*i2)/dble(2*n2);
      dble u3 = dble(-1+n3-2*i3)/dble(2*n3);

      const mat3& bvec = crystal->bvec;
      specialk.push_back( u1*bvec[0] + u2*bvec[1] + u3*bvec[2] );
      reduce_to_bzone(specialk[n++]);
    }
   
    // fold by crystal symmetry operations
    fold_specialk(specialk);
}

void specialk_t::fold_specialk(vector<vec3>& specialk) {
   
    int nspecialk = specialk.size();
    vector<bool> removedk(nspecialk, false);
    
    // for every special k-point not removed
    int n = 0; kstar.clear(); 
    star_t star(crystal);
    for(int i = 0; i < nspecialk; i++) {
      if( removedk[i] ) continue;

      // set the k-point as a representative of the star
      star.k = specialk[i]; 
      kstar.push_back(star); n++;

      // rotate k-point by each crystal symmetry
      // identity should be included in this list
      for(int j = 0; j < crystal->symmetry.size(); j++) {
        vec3 rotk = crystal->symmetry[j] * kstar[n-1].k;

        // see if rotated k-point matches an unremoved k-point
        bool found_match = false;
        for(int l = i; l < nspecialk; l++)
          if( !removedk[l] ) {
            vec3 del = rotk - specialk[l];
            if( mag(del) < tol ) 
              { removedk[l] = true; found_match = true; break; }
          }

        // if so, add corresponding symmetry to star
        if( found_match ) kstar[n-1].symm_index.push_back(j);
      }
    }
    int nstar = n;

    // debug check: print stars
    printf("nstar = %d\n",nstar);
    for(int i = 0; i < nstar; i++) {
      printf("%20.15f %20.15f %20.15f %d\n",
          kstar[i].k.x, kstar[i].k.y, kstar[i].k.z,
          (int)kstar[i].symm_index.size() );
    }
}

void specialk_t::reduce_to_bzone(vec3& kpoint) {

    // for a subset of neighboring reciprocal vectors
    for(int ix = -2; ix <= 2; ix++)
    for(int iy = -2; iy <= 2; iy++)
    for(int iz = -2; iz <= 2; iz++) {
      const mat3& bvec = crystal->bvec;
      vec3 kvec = ix*bvec[0]+iy*bvec[1]+iz*bvec[2];
      if( ix == 0 && iy == 0 && iz == 0 ) continue;

      // find fraction (in r.u.) of k-point outside half-plane 
      dble kfrac;
      kfrac = dot(kpoint-0.5*kvec,kvec);
      kfrac = kfrac/dot(kvec,kvec); 
 
      // perform reduction
      kpoint = kpoint - ceil(kfrac)*kvec;
    }
}


