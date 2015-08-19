/* written by Suffian Khan Aug 2014 */

#include "green0.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <ctime>
#include <map>

#define USE_MATH_DEFINES
#include <cmath>

#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sort.h"
#include "../ext/faddeeva.hh"

using namespace std;
using namespace Faddeeva;

// numerical constants
const double tpi = 2.0*M_PI;
const double fpi = 4.0*M_PI;
const cplx im(0.0,1.0);


// calculate spherical harmonic normalization
// this significantly speeds up Ylm calculation
// note: l increments before m
void green0_t::calc_ylm_norm() {

  // normalization factor sqrt( (2l+1)/4pi (l-m)!/(l+m)! )
  Alm.resize(numLp);
  double fac1 = 1.0;
  for(int m = 0, i = 0; m <= maxlp; m++) {
    double fac2 = 1.0/fac1;
    for(int l = m; l <= maxlp; l++, i++) {
      Alm[i] = sqrt((l+l+1)/(4.0*M_PI));
      Alm[i] *= sqrt(fac2);
      fac2 *= (l+1.0-m)/(l+1.+m);
    }
    fac1 *= (2*m+2)*(2*m+1);
  }
}

// calculate associated Legendre polynomials at cos(th)
// for simplicity and speed, l increments first
void green0_t::calc_plm(double *Plm, 
  const double& cth, const double& sth) {

  // note: avoid integer to double conversion
  double pmm = 1.0, tm1 = 1.0;
  for(int m = 0, i = 0; m <= maxlp; m++) {

    // P[l,l] = (-1)^l (2l-1)!! sin(th)^l
    // P[l+1,l] = (2l+1) cos(th) P[l,l]
    Plm[i] = pmm;
    if(m+1 <= maxlp)
      Plm[i+1] = tm1*cth*pmm;
   
    // (l-m+1)P[l+1,m] = (2l+1) cos(th) P[l,m] - (l+m) P[l-1,m]
    double tl1 = tm1 + 2.0, lm = tm1, lm1 = 2.0;
    for(int l = m+1; l < maxlp; l++, i++) {
      Plm[i+2] = tl1*cth*Plm[i+1]-lm*Plm[i];
      Plm[i+2] /= lm1;
      lm++, lm1++, tl1+=2.0;
    }
    pmm *= -tm1*sth; i+=2;  
    tm1 += 2.0;
  }
}

// calculate v^l Y_lm(v) with normalization from QM
// specifically for use in Ewald summation
// original fast design by William Shelton, & Andrei Smirnov 
void green0_t::calc_vlylm(cplx *vlylm, const vec3& v) {

  // need magnitude and direction vectors 
  double rh2 = v.x*v.x + v.y*v.y;
  double rho = sqrt(rh2);
  double mgv = sqrt(rh2+v.z*v.z);
  double cth = v.z/mgv;
  double sth = rho/mgv;

  // if zero vector return zero array
  if(mgv < 1e-13) {
    for(int L = 0; L < numLp; L++)
      vlylm[L] = cplx(0.0,0.0);
    return;
  }

  // precompute v^l for l = 0..maxl'
  double vl[maxlp+1];
  vl[0] = 1.0; vl[1] = mgv;
  for(int l = 2; l <= maxlp; l++) 
    vl[l] = vl[l-1]*vl[1];

  // generate associated Legendre polynomials
  // for simplicity and speed, l increments first
  double Plm[numLp];
  calc_plm(Plm,cth,sth);

  // explicit representation of exp(i*m*phi) for speed
  const double zerotol = 1e-12;
  double cphi, sphi;
  if(rho > zerotol)
    cphi = v.x/rho, sphi = v.y/rho;
  else
    cphi = 1.0, sphi = 0.0;

  double cmphi = 1.0, smphi = 0.0, cmphi0, smphi0; 
    
  // define v^l Ylm(v) for all L'
  for(int m = 0, i = 0; m <= maxlp; m++) {

    for(int l = m; l <= maxlp; l++, i++) {

      int L0 = l*(l+1);
      double fac1 = vl[l]*Alm[i]*Plm[i]; 

      // circumvent complex arithmetic
      double x =  fac1*cmphi, y = -fac1*smphi;
      vlylm[L0-m] = cplx(x,y);
      if(m % 2 == 0)
        vlylm[L0+m] = cplx(x,-y);
      else
        vlylm[L0+m] = cplx(-x,y);
    }

    // use rotation matrix for speed
    cmphi0 = cmphi, smphi0 = smphi;
    cmphi = cphi*cmphi0 - sphi*smphi0;
    smphi = sphi*cmphi0 + cphi*smphi0;
  }
}

// calculate spherical Hankel function for complex argument
// explicitly for use in real-space structure constants
// assume calculation is up to max l'
void green0_t::calc_hankel(cplx *hl, const cplx& x) {

  // warning: not accurate for |x| << 1

  cplx ix = 1.0/x; 
  cplx ef = exp(im*x)*ix;
  hl[0] = -im*ef;
  hl[1] = -(cplx(1.0,0.0)+im*ix)*ef;
  double tlm1 = 3.0;
  for(int l = 2; l <= maxlp; l++, tlm1 += 2.0)
    hl[l] = tlm1*ix*hl[l-1] - hl[l-2];
}

// perform I(l) = int( x^2l exp(-x^2 r^2-b^2/x^2) )
//   over x = a .. inf 
// a = sqrt(eta)/2, b = i sqrt(z)/2

// Use recursion relation
//   (2l+1) I(l) = 2r^2 I(l+1) - 2b^2 I(l-1) - c(l)
// where
//   c(l)  = a^(2l+1) exp(-a^2 r^2 - b^2/a^2)
//   I(0)  = sqrt(pi)/(4r) [exp(-2br)erfc(ar-b/a)+exp(2br)erfc(ar+b/a)]
//   I(-1) = sqrt(pi)/(4b) [exp(-2br)erfc(ar-b/a)-exp(2br)erfc(ar+b/a)]
// Taken from Abramowitz & Stegun
// Original design by Yang Wang
void green0_t::calc_intfac(cplx *intfac, 
  const cplx& z, const double& r) {

  // warning: need to consider case |b| << 1
  // warning: need to consider case |r| << 1
  // or at least print warning

  // for |r| << 1, simply zero array
  if(r < 1e-13) {
    for(int l = 0; l <= maxlp; l++)
      intfac[l] = cplx(0.0,0.0);
    return;
  }

  // define a and b
  double a = 0.5*sqrt(eta);
  cplx  b = 0.5*im*sqrt(z);

  // setup I(-1) and I(0)
  cplx al = exp(2.0*b*r);
  cplx be = erfc(a*r+b/a);
  cplx ga = erfc(a*r-b/a);
  cplx de = sqrt(M_PI)/4.0;

  cplx Im0 = de/r * (ga/al + al*be);
  cplx Im1 = de/b * (ga/al - al*be);

  // save I(0) and calculate I(1)
  al = 0.5/(r*r), be = 2.0*b*b, ga = a*a;
  de = exp(-ga*r*r-b*b/ga);
  cplx c = a*de; 
  intfac[0] = Im0;
  if(maxlp == 0) return;
  intfac[1] = al*(Im0 + be*Im1 + c);

  // use recursion relation to fill in l > 0
  c *= ga;
  for(int l = 2; l <= maxlp; l++) { 
    intfac[l] = (2.*l-1.)*intfac[l-1] + be*intfac[l-2] + c;
    intfac[l] *= al; c *= ga;
  }
}


