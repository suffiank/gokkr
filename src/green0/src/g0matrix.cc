#include "green0.h"
#include <cstdio>
using namespace std;

void green0_t::calc_g0mat(cplx *Gmat, const cplx& z, const vec3& k) {

  // note: it was more convenient to have one routine that accesses
  // k-vectors on an indexed list than separate routines for both cases
 
  // add 'k' to the list of kpoints
  kpoints.push_back(k);
  int kindex = kpoints.size()-1;

  // if caching enabled, precompute (k+K)^l Y_lm(k+K)
  // note: the below operations do not conserve memory
  if(cachemode & CACHE_KSPACE) {
    vector<cplx> ylm(numLp);
    vector< vector<cplx> > klylm(klatt.size(),ylm);
    klylm_latt.push_back(klylm);
    for(int i = 0; i < klatt.size(); i++) {
      vec3 kvec = k + klatt[i];
      calc_vlylm(klylm_latt[kindex][i].data(),kvec);
    }
  }

  // call back-end for G0 calculation with new kpoint
  calc_g0mat(Gmat, z, kindex);
  
  // remove kpoint (and cache) from list to prevent side-effects
  kpoints.pop_back();
  if(cachemode & CACHE_KSPACE) klylm_latt.pop_back();
}

void green0_t::calc_g0mat(cplx* Gmat, const cplx& z, const int kindex) {

  int natom = basis.size();
  int matdim = natom*numL;

  // begin timer if clocking enabled
  if(log_times) timer.begin("Generate G0 matrix");

  // for every unique aij vector
  cplx Gij[numL*numL], Gji[numL*numL];
  for(int n = 0; n < aijlist.size(); n++) {

    // determine source (i0,j0) for this aij set
    int i0 = aij_sources[n][0].first;
    int j0 = aij_sources[n][0].second;

    // fill-in structure constants for aij
    cplx *Gmat_i0j0 = &Gmat[0] + i0*numL*matdim + j0*numL;
    calc_g0ij_kspace(Gmat_i0j0, matdim, z, n, kindex);

    // for every other corresponding (i,j) atom pair
    for(int m = 1; m < aij_sources[n].size(); m++) {

      // determine (i,j) for this pair
      int i = aij_sources[n][m].first;
      int j = aij_sources[n][m].second;
 
      // copy to appropriate position in larger matrix
      cplx *Gmat_ij = &Gmat[0] + i*numL*matdim + j*numL;
      for(int L1 = 0; L1 < numL*matdim; L1+=matdim)
      for(int L2 = 0; L2 < numL; L2++)
        Gmat_ij[L1+L2] = Gmat_i0j0[L1+L2];
    }
  }

  // print time if clocking enabled
  if(log_times) timer.end();
}


