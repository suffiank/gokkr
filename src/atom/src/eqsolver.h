#include <cstdio>
#include <cstdlib>

#include "../../util/src/mathfn.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_sf_bessel.h"
using namespace std;

// Find solution of
//   [-(d/dr)^2 + l(l+1)/r^2 + V(r)]P_nl(r) = E_nl P_nl(r)
// where P_nl(r) = r R_nl(r)
 
// Require change of independent variable to x = log(r/mu):
//   [-(d/dx)^2 + d/dx + l(l+1) + r^2 (V-E)] P = 0

// Define Q = dP/dx to change to a first order equation:
//   d[P, Q]/dx = [Q, fP+Q] 
// where f(x) = l(l+1) + r(x)^2 ( V(x)-E )


// define taylor series expansion to O(r^5) about r = 0: 
//   P(x) = sum(n=0..4, c_n r^n)
// warning: may not be working to correct order
template<typename T>
int solve_sch_taylor_origin(const xfunc& r, func<T>& P, func<T>& Q, 
    int l, T f1, T f2, interval_t range) {

  // set boundary condition P_nl(r << 1) = A r^(l+1)
  // recall index=0 corresponds to r=r_min, and Q = dP/dx 
  // note: can we generate an improved expansion based on potential?
  T c[5];
  c[0] = 1.0; 
  c[1] = f1/(2.*(l+1.0)); 
  c[2] = (c[1]*f1+f2)/(2.*(l+3.0));
  
	// fill in P and Q for four indices using taylor
  // note: can multiplying by r^(l+1) lead to round-off?
  for(int i = range.a; i <= range.b; i++) {
		P[i] = 0.0; Q[i] = 0.0;
    for(int k = 2; k >= 0; k--) {
			P[i] =   c[k] + r[i]*P[i];
			Q[i] = (k+l+1.)*c[k] + r[i]*Q[i];
		}

    for(int k = 0; k < l+1; k++)
      { P[i] *= r[i]; Q[i] *= r[i]; }
  }

}

// fill in tail using WKB approximation
template<typename T>
void solve_sch_wkb_tail(dble dx, const xfunc& r, func<T>& P, func<T>& Q, 
    const func<T>& g, interval_t range) {

  // note: routine only valid for exponential decay 
  int WKB = range.a;
  T A = 1.0e-03*sqrt(sqrt(g[WKB]));
  T kappa = 0.0;
  for(int i = range.a; i <= range.b; i++) {
  
     dble ep = der(dx, g, i); 
     P[i] = A*exp(-kappa)/sqrt(sqrt(g[i]));
     Q[i] = -(sqrt(g[i])*r[i]+ep/(4.0*g[i]))*P[i];
    
     if( i < range.b ) // trapezoidal rule (do better?) 
       kappa += 0.5*(sqrt(g[i+1])+sqrt(g[i]))*(r[i+1]-r[i]);
  }
}

template<int step = 1, typename T>
int solve_sch_rk4(dble dx, func<T>& P, func<T>& Q, const func<T>& f, interval_t range) {

  // use template so 'step' is a compile-time constant
  // use step = +1 for forward solve
  // use step = -1 for backward solve
  if( step != 1 && step != -1 )
    { std::printf("solve_sch_abm4:bad step"); exit(1); }

  // use GSL to interpolate f(x) for midpoints
  // todo: slow! compute this only once for all n,l?
  /* xfunc x(f.size());
  for(int i = 0; i < f.size(); i++) x[i] = i*dx;
  gsl_interp *pfi = gsl_interp_alloc(gsl_interp_cspline,f.size());
  gsl_interp_accel *pfa = gsl_interp_accel_alloc();
  gsl_interp_init(pfi, x.data(), f.data(), f.size()); */

  T k1P,k1Q,k2P,k2Q,k3P,k3Q,k4P,k4Q,f23;
 
  int nodes = 0;
  dble sign = 1.0;

  dx = step*dx;
  for(int i = range.a; i != range.b+step; i += step) {
 
    // get half-point using some interpolation method
    f23 = (f[i-step]+f[i])/2.0;
    // f23 = gsl_interp_eval(pfi, x.data(), f.data(), (i-0.5*step)*(step*dx), pfa);

    // advance point using Runge-Kutta scheme
    k1P = dx*Q[i-step];           k1Q = dx*f[i-step]*P[i-step] + k1P;
    k2P = dx*(Q[i-step]+k1Q/2.0); k2Q = dx*f23*(P[i-step]+k1P/2.0) + k2P;
    k3P = dx*(Q[i-step]+k2Q/2.0); k3Q = dx*f23*(P[i-step]+k2P/2.0) + k3P;
    k4P = dx*(Q[i-step]+k3Q);     k4Q = dx*f[i]*(P[i-step]+k3P) + k4P;
    P[i] = P[i-step] + (k1P + 2.0*k2P + 2.0*k3P + k4P)/6.0;
    Q[i] = Q[i-step] + (k1Q + 2.0*k2Q + 2.0*k3Q + k4Q)/6.0;

    // check for zero crossing
    if( real_part(P[i])*sign < 0.0 )
      { nodes++; sign = -sign;}
  }

  // free GSL resources
  // gsl_interp_free(pfi); gsl_interp_accel_free(pfa);

  return nodes;
}

template<int step = 1, typename T>
int solve_sch_abm4(dble dx, func<T>& P, func<T>& Q, const func<T>& f, interval_t range) {

  if( step != 1 && step != -1 )
    { std::printf("solve_sch_abm4:bad step"); exit(1); }

  // predictor-corrector coefficients
  const dble a[4] =  {55./24.*dx,-59./24.*dx,37./24.*dx,-9./24.*dx};
  const dble b[4] =  {3./8.*dx,19./24.*dx,-5./24.*dx,1./24.*dx};
  const int ncorr = 3;
  T fP[4], dP, dQ;
 
  int nodes = 0;
  dble sign = 1.0;

  // setup memory function for f[j]*P[j]
  const int i0 = range.a;
  const int di = step;
  fP[0] = f[i0-1*di]*P[i0-1*di];
  fP[1] = f[i0-2*di]*P[i0-2*di];
  fP[2] = f[i0-3*di]*P[i0-3*di];
  fP[3] = f[i0-4*di]*P[i0-4*di];

  // proceed with 4-step Adams-Bashforth-Moulton integrator
  for(int i = range.a; i != range.b+step; i += step) {

    // predictor
    dP = a[0]*Q[i-di] + a[1]*Q[i-2*di] + a[2]*Q[i-3*di] + a[3]*Q[i-4*di];
    dQ = a[0]*fP[0] + a[1]*fP[1] + a[2]*fP[2] + a[3]*fP[3] + dP;
    P[i] = P[i-di] + double(step)*dP; Q[i] = Q[i-di] + double(step)*dQ;
    
    // corrector
    for(int k = 0; k < ncorr; k++) {
      dP = b[0]*Q[i] + b[1]*Q[i-di] + b[2]*Q[i-2*di] + b[3]*Q[i-3*di];
      dQ = b[0]*f[i]*P[i] + b[1]*fP[0] + b[2]*fP[1] + b[3]*fP[2] + dP;
      P[i] = P[i-di] + double(step)*dP; Q[i] = Q[i-di] + double(step)*dQ; 
    }

    // update array f P to coincide with last four values
    fP[3] = fP[2]; fP[2] = fP[1]; fP[1] = fP[0]; fP[0] = f[i]*P[i];

    // check for zero crossing
    if( real_part(P[i])*sign < 0.0 )
      { nodes++; sign = -sign;}
  }

  return nodes;
}

// fill in ricati-bessel r f_l(pr) 
enum bessel_t { bessel_j, bessel_y, bessel_h1, bessel_h2 };
template<bessel_t filltype, typename T>
void solve_sch_riccati(T E, int l, const xfunc& r, func<T>& P, func<T>& Q,
    interval_t range) {

  // fl0 = f_l(pr), fl1 = f_(l+1)(pr), flp = f_l'(pr)
  T p = sqrt(E);
  for(int i = range.a; i <= range.b; i++) {

    T fl0, fl1;
    switch(filltype) {
      case bessel_j: 
        fl0 = sph_bessel_jl(l,   p*r[i]);
        fl1 = sph_bessel_jl(l+1, p*r[i]);
        break;
      case bessel_y: 
        fl0 = sph_bessel_yl(l,   p*r[i]);
        fl1 = sph_bessel_yl(l+1, p*r[i]);
        break;
      case bessel_h1: 
        fl0 = sph_bessel_hl1(l,   p*r[i]);
        fl1 = sph_bessel_hl1(l+1, p*r[i]);
        break;
      default:
        throw "solve_sch_sph_bessel :: invalid bessel type";
    }

    // f_l'(z) = -f_(l+1)(z) + l/z f_l(z)
    // see NIST 'digital library of mathematical functions'
    T flp = -fl1 + dble(l)/(p*r[i])*fl0;

    P[i] = r[i]*fl0;
    Q[i] = r[i]*fl0 + p*r[i]*r[i]*flp;
  }
}
