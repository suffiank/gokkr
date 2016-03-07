#ifndef _MATHUTIL_H
#define _MATHUTIL_H

#include <cmath>
#include <vector>
#include <utility>
#include <complex>

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))

// type definitions
typedef double dble;
typedef std::complex<dble> cplx;

inline dble real_part(dble x) {return x;}
inline dble real_part(cplx z) {return z.real();}
inline dble imag_part(cplx z) {return z.imag();}

template<typename T>
using func = std::vector<T>;
typedef func<dble> xfunc;
typedef func<cplx> zfunc;

struct interval_t { 
  int a,b; 
  interval_t(int aa, int bb) : a(aa), b(bb) {}
};

// constants 
const dble pi = M_PI;
const dble tpi = 2.0*M_PI;
const dble fpi = 4.0*M_PI;
const cplx im = cplx(0.0,1.0);

// fitting
// warning: poor precision unless near x = 0
template<typename T, typename F>
void getcubicfit(const F& x, const F& f, T *a)
{
  T al,be,ga,de,ep;

  al = f[0]/( (x[0]-x[1])*(x[0]-x[2])*(x[0]-x[3]) );
  be = f[1]/( (x[1]-x[0])*(x[1]-x[2])*(x[1]-x[3]) );
  ga = f[2]/( (x[2]-x[0])*(x[2]-x[1])*(x[2]-x[3]) );
  de = f[3]/( (x[3]-x[0])*(x[3]-x[1])*(x[3]-x[2]) );

  a[3] = al + be + ga + de;

  ep = x[0]+x[1]+x[2]+x[3];
  a[2] = al*( x[0] - ep ) + 
         be*( x[1] - ep ) + 
         ga*( x[2] - ep ) + 
         de*( x[3] - ep );

  a[1] = al*(x[1]*x[2]+x[1]*x[3]+x[2]*x[3]) + 
         be*(x[0]*x[2]+x[0]*x[3]+x[2]*x[3]) + 
         ga*(x[0]*x[1]+x[0]*x[3]+x[1]*x[3]) + 
         de*(x[0]*x[1]+x[0]*x[2]+x[1]*x[2]);

  a[0] = -al*( x[1]*x[2]*x[3] ) 
         -be*( x[0]*x[2]*x[3] ) 
         -ga*( x[0]*x[1]*x[3] ) 
         -de*( x[0]*x[1]*x[2] );
}

// derivatives (first, second, left, right)

// 7-point central difference
template<typename T>
T der(dble dx, const func<T>& f, int i) {
 return (-1./60.*f[i-3]+3./20.*f[i-2]-3./4.*f[i-1]+0.*f[i] \
   +3./4.*f[i+1]-3./20.*f[i+2]+1./60.*f[i+3])/dx;
}

template<typename T>
T der2(dble dx, const func<T>& f, int i) {
 return (1./90.*f[i-3]-3./20.*f[i-2]+3./2.*f[i-1]-49./18.*f[i] \
   +3./2.*f[i+1]-3./20.*f[i+2]+1./90.*f[i+3])/(dx*dx); 
}

// 7-point left derivative
template<typename T>
T lder(dble dx, const func<T>& f, int i) {
  return (49./20.*f[i]-6.*f[i-1]+15./2.*f[i-2]-20./3.*f[i-3] \
    +15./4.*f[i-4]-6./5.*f[i-5]+1.0/6.*f[i-6])/dx;
}

template<typename T>
T lder2(dble dx, const func<T>& f, int i) {
  return (469./90.*f[i]-223./10.*f[i-1]+879./20.*f[i-2] \
    -949./18.*f[i-3]+41.*f[i-4]-201./10.*f[i-5] \
    +1019./180.*f[i-6]-7./10.*f[i-7])/(dx*dx);
}

// 7-point right derivative
template<typename T>
T rder(dble dx, const func<T>& f, int i) {
  return (-49./20.*f[i]+6.*f[i+1]-15./2.*f[i+2]+20./3.*f[i+3] \
    -15./4.*f[i+4]+6./5.*f[i+5]-1.0/6.*f[i+6])/dx;
}

template<typename T>
T rder2(dble dx, const func<T>& f, int i) {
  return (469./90.*f[i]-223./10.*f[i+1]+879./20.*f[i+2] \
    -949./18.*f[i+3]+41.*f[i+4]-201./10.*f[i+5] \
    +1019./180.*f[i+6]-7./10.*f[i+7])/(dx*dx);
}

// integration routines 
template<typename T>
void trapint(dble dx, const func<T>& f, func<T> *I) {
  (*I)[0] = 0.0;
  for(int i = 1; i < f.size(); i++)
    (*I)[i] = (*I)[i-1] + 0.5*(f[i]+f[i-1])*dx;
}

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

// root-mean-square of two functions
// void getrms(const func& f_in, const func& f_out);

// spherical bessel functions
template<typename T>
T sph_bessel_jl(int l, T x) {

  T j[l+1];

  // use series expansion for |x| << 1 
  // jl(x) = x^l sum( (-1/2 x^2)^k / k!(2l+2k+1)!!, k = 0,+inf )
  if( abs(x) < 0.001 ) {

    int nterm = 3;

    dble f1 = 1.0, f2 = 1.0;
    for(int i = 2*l+1; i > 0; i-=2) 
      f2 *= i;

    j[l] = 1.0;
    for(int i = 0; i < nterm; i++) {
      j[l] += pow(-0.5*x*x,i)/f1/f2; 
      f1 *= (i+1.); f2 *= (2.*l+2.*i+3.);
    }
    j[l] *= pow(x,l);
    return j[l];
  }

  // recursion relations for |x| != 0
  // j_(l-1) + j_(l+1) = (2l+1) j_l/x 
  j[0] = sin(x)/x;
  if( l > 0 ) j[1] = sin(x)/(x*x) - cos(x)/x;

  for(int i = 2; i <= l; i++)
    j[i] = (2.*i-1.)*j[i-1]/x - j[i-2];

  return j[l];
}

template<typename T>
T sph_bessel_yl(int l, T x) {

  T y[l+1];
  y[0] = -cos(x)/x;
  if( l > 0 ) y[1] = -cos(x)/(x*x) - sin(x)/x;

  for(int i = 2; i <= l; i++)
    y[i] = (2.*i-1.)*y[i-1]/x - y[i-2];

  return y[l];
}

template<typename T>
T sph_bessel_hl1(int l, T x) {

  T h[l+1];
  h[0] = -im*exp(im*x)/x;
  if( l > 0 ) h[1] = -(x+im)/(x*x)*exp(im*x);

  for(int i = 2; i <= l; i++)
    h[i] = (2.*i-1.)*h[i-1]/x - h[i-2];

  return h[l];
}

#endif
