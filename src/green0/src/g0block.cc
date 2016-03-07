/* written by Suffian Khan Aug 2014 */

#include "green0.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <ctime>
#include <map>

#define USE_MATH_DEFINES
#include <cmath>
#include "../../util/src/ylm.h"

using namespace std;

// calculation of K-SPACE structure constants 'B^ij_LiLj' 
// see Electron Scattering by Zabloudil, pg. 158 
// note: expression to match MECCA is different!
void green0_t::calc_g0ij_kspace(cplx* Gij, const int stride, 
  const cplx& z, const int rindex, const int kindex) {

  // relabel 'r' and 'k' vectors for conveniance
  const vec3& r = aijlist[rindex];
  const vec3& k = kpoints[kindex];

  // if energy changed, recalculate integration factors
  if( z != lastz ) {
    if(cachemode & CACHE_RSPACE) 
      if(calcmode == USE_EWALD)
        calc_all_intfac(z); 
      else
        calc_all_hankel(z);
    d00 = calc_d00(z);
    lastz = z;
  }   

  // define el- momentum
  cplx p = sqrt(z);
  if( p.imag() < 0.0 ) p = -p;

  // zero Dlm array
  cplx Dlm[numLp];
  for(int L = 0; L < numLp; L++) 
    Dlm[L] = cplx(0.0,0.0); 
 
  if(calcmode == USE_EWALD) { 

    // perform dlm calculation via Ewald method

    // perform r-space lattice summation
    if(cachemode & CACHE_RSPACE)
      add_dlm_ewald_rsum(Dlm, z, p, rindex, k);
    else
      add_dlm_ewald_rsum(Dlm, z, p, r, k);   
  
    // perform k-space lattice summation
    if(cachemode & CACHE_KSPACE) 
      add_dlm_ewald_ksum(Dlm, z, p, r, kindex);
    else
      add_dlm_ewald_ksum(Dlm, z, p, r, k);
  
    // add D^(3)_00 if on diagonal block
    if( mag(aijlist[rindex]) <= 1e-13 )
      Dlm[0] += d00;
  }
  else { 

    // perform dlm calculation via Fourier transform
    if(cachemode & CACHE_RSPACE)
      calc_dlm_fourier(Dlm, z, p, rindex, k);
    else
      calc_dlm_fourier(Dlm, z, p, r, k);
  }
  
  // perform gaunt summation
  for(int L1 = 0, L12 = 0; L1 < numL*stride; L1 += stride)
  for(int L2 = 0; L2 < numL; L2++, L12++) {
    Gij[L1+L2] = cplx(0.0,0.0);
    for(int i = gindex[L12]; i != gindex[L12+1]; i++)
      Gij[L1+L2] += gaunt[i]*Dlm[gLpval[i]];
  }

  // final factor of 4pi i^(l-l') exp(-ik.r) (-1)^(m+m')
  // note: minus signs introduced to agree with MECCA?
  cplx fac1 = fpi*exp(cplx(0.0,-dot(k,r))), fac2;
  for(int l1 = 0, L1 = 0; l1 <= maxl; l1++, fac1 *= im)
  for(int m1 = -l1; m1 <= l1; m1++, L1 += stride, fac1 = -fac1) {
    fac2 = 1.0;
    for(int l2 = 0, L2 = 0; l2 <= maxl; l2++, fac2 /= im) 
    for(int m2 = -l2; m2 <= l2; m2++, L2++, fac1 = -fac1)
      Gij[L1+L2] *= fac1*fac2;
  }
} 


// perform sqrt(eta)/2pi sum( (z/eta)^n/[(2n-1)n!], n=0..inf )
inline cplx green0_t::calc_d00(const cplx& z) {

  // warning: may suffer inaccuracy for high terms?

  const double tol = 1.e-14;
  const int nmin = 10;

  cplx zoe = z/eta;
  cplx zoen = 1.0;
  
  int n = 0; double err = 1.e10;
  cplx sum = 0.0, nf = cplx(1.0,0.0);
  while(n <= nmin || err > tol) {
    cplx term = zoen/( (2.*n-1)*nf );
    sum += term; 
    err = abs(term/sum);
    n++; nf = cplx(n,0.0)*nf; zoen *= zoe;
  }
  sum *= -sqrt(eta)/(2.*M_PI);
  return sum;
}

inline void green0_t::add_dlm_ewald_rsum(cplx *Dlm, const cplx& z, 
  const cplx& p, const int rindex, const vec3& k) {

  // perform r-space lattice summation:
  // -2/sqrt(pi) (-2i/p)^l *
  //   sum( all R, |R+r|^l Y_lm(R+r) exp(ik.(R+r) ) *
  //   int(x^2l exp(-x^2 (R+r)^2 + z/4x^2)dx ) )

  const vec3& r = aijlist[rindex];

  cplx A = -2.0/sqrt(M_PI)*exp(cplx(0.0,dot(k,r))); 
  cplx c = -cplx(0.0,2.0)/p;
  for(int i = 0; i < rlatt.size(); i++) {

    cplx fac = A*exp(cplx(0.0,dot(rlatt[i],k)));
    for(int l = 0, L = 0; l <= maxlp; l++) {
      for(int m = -l; m <= l; m++, L++)
        Dlm[L] += fac*intfac[rindex][i][l]*rlylm_latt[rindex][i][L];
      fac *= c;
    }
  }

}

// same as above except no cached harmonics
inline void green0_t::add_dlm_ewald_rsum(cplx *Dlm, const cplx& z, 
  const cplx& p, const vec3& r, const vec3& k) {

  // todo: do one shell at a time
  cplx rlylm[numLp], intfac[maxlp+1];
  cplx A = -2.0/sqrt(M_PI)*exp(cplx(0.0,dot(k,r))); 
  cplx c = -cplx(0.0,2.0)/p;
  for(int i = 0; i < rlatt.size(); i++) {

    vec3 rv = rlatt[i] + r;
    const double mrv = mag(rv);
    calc_vlylm(maxlp, Alm.data(), rlylm,rv);
    calc_intfac(intfac,z,mrv);

    cplx fac = A*exp(cplx(0.0,dot(rlatt[i],k)));
    for(int l = 0, L = 0; l <= maxlp; l++) {
      for(int m = -l; m <= l; m++, L++)
        Dlm[L] += fac*intfac[l]*rlylm[L];
      fac *= c;
    }
  }
}

inline void green0_t::add_dlm_ewald_ksum(cplx *Dlm, const cplx& z, 
  const cplx& p, const vec3& r, const vec3& k) {

  // perform k-space lattice summation:
  // -(4pi/volume) exp(z/eta) p^(-l) *
  //   sum( all K, exp(-(K+k)^2/eta)/(|K+k|^2-z) *
  //    |K+k|^l Y_lm(K+k) exp(iK.r) )

  // todo: do one shell at a time
  cplx klylm[numLp];
  cplx A = -fpi/vol*exp(z/eta), c = 1.0/p;
  for(int i = 0; i < klatt.size(); i++) {

    vec3 kv = klatt[i] + k;
    double kv2 = dot(kv,kv);
    calc_vlylm(maxlp,Alm.data(),klylm,kv);

    cplx fac = A;
    fac *= exp(cplx(-kv2/eta,-dot(klatt[i],r)));
    fac /= kv2-z; 

    for(int l = 0, L = 0; l <= maxlp; l++) {
      for(int m = -l; m <= l; m++, L++) 
        Dlm[L] += fac*klylm[L];
      fac *= c;
    }
  } 
}

// same as above except |K+k|^l Y_lm(K+k) is cached
inline void green0_t::add_dlm_ewald_ksum(cplx *Dlm, const cplx& z, 
  const cplx& p, const vec3& r, int kindex) {

  //  original version
  cplx A = -fpi/vol*exp(z/eta), c = 1.0/p;
  for(int i = 0; i < klatt.size(); i++) {

    vec3 kv = klatt[i] + kpoints[kindex];
    double kv2 = dot(kv,kv);

    cplx fac = A;
    fac *= exp(cplx(-kv2/eta,-dot(klatt[i],r)));
    fac /= kv2-z; 

    for(int l = 0, L = 0; l <= maxlp; l++) {
      for(int m = -l; m <= l; m++, L++) 
        Dlm[L] += fac*klylm_latt[kindex][i][L];
      fac *= c;
    }
  } 

}

inline void green0_t::calc_dlm_fourier(cplx *Dlm, const cplx& z, 
  const cplx& p, const vec3& r, const vec3& k) {

  for(int L = 0; L < numLp; L++)
    Dlm[L] = 0.0;

  if(mag(r) < 1.e-10) 
    Dlm[0] = -im*p/sqrt(fpi);

  for(int i = 0; i < rlatt.size(); i++) {
 
    vec3 rvec = r + rlatt[i];
    double mgr = mag(rvec);
    if(mgr < 1.e-10) continue;

    // get h_l(pr) Y_lm(r) for spherical Hankel h_l(x)
    // note: not possible to use GSL due to complex argument
    cplx hl[maxlp+1];
    calc_hankel(hl, p*mag(rvec));
  
    cplx ylm[numLp];
    calc_vlylm(maxlp, Alm.data(), ylm, rvec);
    double mrv = mag(rvec);
  
    cplx fac = -im, expf = p*exp(cplx(0.0,dot(k,rvec)));
    for(int l = 0, L = 0; l <= maxlp; l++, fac *= -im/mrv)
    for(int m = -l; m <= l; m++, L++)
      Dlm[L] += fac*expf*hl[l]*ylm[L];
  } 
}

// overloaded, cached version of above
inline void green0_t::calc_dlm_fourier(cplx *Dlm, const cplx& z, 
  const cplx& p, const int rindex, const vec3& k) {

  for(int L = 0; L < numLp; L++)
    Dlm[L] = 0.0;

  if(rindex == 0) 
    Dlm[0] = -im*p/sqrt(fpi);

  for(int i = 0; i < rlatt.size(); i++) {
 
    vec3 rvec = aijlist[rindex] + rlatt[i]; 
    double mrv = mag(rvec);
    if(mrv < 1.e-10) continue;

    // note: a factor of 1/(r+R)^l to cancel that in (r+R)^l Y_lm(r+R)
    cplx fac = -im, expf = p*exp(cplx(0.0,dot(k,rvec)));
    for(int l = 0, L = 0; l <= maxlp; l++, fac *= -im/mrv )
    for(int m = -l; m <= l; m++, L++)
      Dlm[L] += fac*expf*hankel[rindex][i][l]*rlylm_latt[rindex][i][L];
  } 
}


