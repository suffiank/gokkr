#include "mathfn.h"

// Gauss-Legendre table
void fill_gauss_legendre_table(int n, dble* x, dble* w) {

  const double tol = 1.e-15;

  // see Numerical Recipes
  // generate ith Gauss-Legendre root in +x and -x pairs
  const int m = n/2 + 1;
  for(int i = 0; i < m; i++) {

    // starting guess for Legendre root
    dble z = cos(pi*(i+0.75)/(n+0.5));

    recalculate:

    // recursive calculation of Legendre
    //   (n+1) P[n+1](x) = (2n+1) x P[n](x) - n P[n-1](x)
    dble k = 0.0, Pk = 1.0, Pk1 = 0.0, Pk2;
    for(int j = 0; j < n; j++) {
      Pk2 = Pk1; Pk1 = Pk; k++;
      Pk = ((2.0*k-1.0)*z*Pk1 - (k-1.0)*Pk2)/k;
    } 
    
    // Newton's method to improve root
    //   x' = x - P'[n](x)/P[n](x)
    // where Legendre P'[n](x) satisfies
    //   (x^2 - 1)/n P[n]'(x) = x P[n](x) - P[n-1](x)
    dble dPk = (z*Pk - Pk1)*k/(z*z-1.0);
    dble zp = z - Pk/dPk;
    if( abs(zp-z) > tol )
      { z = zp; goto recalculate; }

    // save ith Gauss-Legendre root and weight
    //  weight = 2/[ (1-x^2) P'[n](x)^2 ]
    x[i] = z; w[i] = 2./( (1.-z*z)*dPk*dPk );
    x[n-1-i] = -x[i]; w[n-1-i] = w[i];
  }
}


