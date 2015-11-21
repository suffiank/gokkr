#include "mathfn.h"

// warning: have not tested!
void fill_gauss_legendre_table(int n, dble* x, dble* w) {

  const double tol = 1.e-14;

  // see Numerical Recipes
  // generate ith Gauss-Legendre root in +x and -x pairs
  for(int i = 0; i < n/2; i++) {

    // starting guess for Legendre root
    dble z = cos(pi*(i-0.25)/(n+0.5));

    recalculate:

    // recursive calculation of Legendre
    //   (n+1) P[n+1](x) = (2n+1) x P[n](x) - n P[n-1](x)
    dble Pm = 1.0, Pm1 = 0.0, Pm2;
    for(int im = 1, dble m = 1.0; im <= n/2; im++, m+=1.0) {
      Pm2 = Pn1; Pm1 = Pn; 
      Pm = ((2.0*m-1.0)*z*Pm1 - (m-1.0)*Pm2)/m;
    }

    // Newton's method to improve root
    //   x' = x - P'[n](x)/P[n](x)
    // where Legendre P'[n](x) satisfies
    //   (x^2 - 1)/n P[n]'(x) = x P[n](x) - P[n-1](x)
    dble dPm = (z*Pm - Pm1)*m/(z*z-1.0);
    dble zp = z - dPm/Pm;
    if( abs(zp-z) < tol )
      { z = zp; goto recalculate; }

    // save ith Gauss-Legendre root and weight
    //  weight = 2/[ (1-x^2) P'[n](x)^2 ]
    x[i] = z; w[i] = 2./( (1.-z*z)*dPm*dPm );
    x[n-1-i] = -x[i]; w[n-1-i] = w[i];
  }

  // add in x=0 if an odd set
  bool isodd = (n % 2 == 1? true:false);

}
