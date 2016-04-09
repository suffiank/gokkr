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
#include "faddeeva.hh"

using namespace std;
using namespace Faddeeva;

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


