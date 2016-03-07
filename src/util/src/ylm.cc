#include "ylm.h"

// calculate spherical harmonic normalization
// Alm = sqrt( (2l+1)/4pi (l-m)!/(l+m)! ) 
// this significantly speeds up Ylm calculation
// note: l increments before m
void calc_ylm_norm(const int maxl, double *Alm) {

  const int& maxlp = maxl;

  // normalization factor sqrt( (2l+1)/4pi (l-m)!/(l+m)! )
  double fac1 = 1.0;
  for(int m = 0, i = 0; m <= maxlp; m++) {
    double fac2 = 1.0/fac1;
    for(int l = m; l <= maxlp; l++, i++) {
      Alm[i] = sqrt((l+l+1)/fpi);
      Alm[i] *= sqrt(fac2);
      fac2 *= (l+1.0-m)/(l+1.+m);
    }
    fac1 *= (2*m+2)*(2*m+1);
  }
}

// calculate associated Legendre polynomials at cos(th)
// for simplicity and speed, l increments first
void calc_plm(const int maxl, double *Plm, const double& cth, const double& sth) {

  const int& maxlp = maxl;

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
void calc_vlylm(const int maxl, const double *Alm, cplx *vlylm, const vec3& v) {

  const int& maxlp = maxl;
  const int numLp = (maxlp+1)*(maxlp+1);

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
  calc_plm(maxlp,Plm,cth,sth);

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


