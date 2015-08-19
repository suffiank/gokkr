/* 

A version of the Herman & Skillman program for the non-relativistic
Hartree-Fock-Slater solutions of a single atom 

Started by Suffian Khan, 04/24/2011 
First completion, 05/10/2015

exclusively using atomic Rydberg units
(i.e. hbar = 2m = e^2/2 = a0 = 1, c = 2/alpha)

*/

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

int main(int argc, char *argv[])
{

  // using atomic Rydberg units	
  int   Z = 7;    // atomic number
  int   N;         // number of grid positions
  real* r;         // radial grid positions
  real* rV;        // radius*potential at grid positions
  real* si;        // 4*pi*r*r*(n(r) = el- density) on grid
  real* P;
  real* Q;
  real* f;
  real* g;

  if( argc == 1 )
    N = 10000;
  else
    N = atoi(argv[1]);

  int*   orb_n;  // orbital principle quantum numbers
  int*   orb_l;  // orbital angular momentums
  int*   orb_w;  // orbital el- occupancy
  real*  orb_E;  // orbital energies, in order of filling
  real** orb_S;  // orbital sigma = 4*pi*r*r*n(r)
  real** orb_P;  // orbital solution P_nl(r) = r R_nl(E,r)

  int i,j,k;

  // setup a logarithmic grid
  // use Thomas-Fermi universal scaling factor: 
  //   r = mu r_TF = (1/2)(3pi/4)^(2/3) Z^(-1/3) r_TF

  // const real xi = log(0.001);
  const real xi = log(0.0001);
  const real dx = ( log(200.0)-xi )/(N-1);
  const real xf = xi + (N-1)*dx;  
  const real mu = 0.88534138*pow(Z,-1.0/3.0);

  r = new real[N];
  real x = xi;
  for(i = 0; i < N; i++)
    r[i] = mu*exp(xi + i*dx);

  // fill in Latter's fit of Thomas-Fermi potential
  // see Phys.Rev. v99 p510 by Latter
  // see DFT by Parr & Wang, Sec 6.2 pg 110 (in Hartree units)

  rV = new real[N];
  si = new real[N];
  P = new real[N];
  Q = new real[N];
  f = new real[N];
  g = new real[N];
  real rx;
  for(i = 0; i < N; i++) {
    x = r[i]/mu; rx = sqrt(x); 
    rV[i] = 1.0 + 0.02747*rx + 1.243*x - 0.1486*x*rx + 0.2302*x*x + \
            0.007298*x*x*rx + 0.006944*x*x*x;
    rV[i] = 2.0*Z/rV[i];
    si[i] = 4.0/(3.0*pi)*rV[i]*sqrt(rV[i]*r[i]);
    rV[i] = -rV[i];
  }

  // replace with standard potential for testing energy levels
  // for(i = 0; i < N; i++)
  //  rV[i] = -2.0*Z; 
   
  // replace with spherical square well 
  real U0 = -10.0, R0 = 1.0;
  for(i = 0; i < N; i++)
    rV[i] = (r[i] < R0? r[i]*U0 : 0.0); 

  // replace with isotropic harmonic oscillator
  /* real omega = 1.0;
  for(i = 0; i < N; i++)
    rV[i] = 0.25*omega*omega*r[i]*r[i]*r[i]; */
  
  // use Madelung's rule to guess highest filled nl subshell

  int n, l, nlvl = 0, nfill = 0;
  for(i = 1; i <= 8; i++)
  for(n = i/2+1, l = (i-1)/2; n <= i; n++, l--)
  if( nfill < Z ) 
    { nfill += MIN(2*(2*l+1), Z-nfill); nlvl++; }
 
  // fill in a guess electronic configuration
  // keep one extra level in case filling changes
  // take initial energy guess as from V(r) = -2Z/r

  nlvl++;
  orb_n = new int[nlvl];
  orb_l = new int[nlvl];
  orb_w = new int[nlvl];
  orb_E = new real[nlvl];
  orb_S = new real*[nlvl];
  orb_P = new real*[nlvl];
  for(i = 0; i < nlvl; i++) {
    orb_S[i] = new real[N];
    orb_P[i] = new real[N];
  }

  j = 0, nfill = 0;
  for(i = 1; i <= 8; i++)
  for(n = i/2+1, l = (i-1)/2; n <= i; n++, l--)
  if( j < nlvl ) {
    orb_n[j] = n; orb_l[j] = l;
    orb_E[j] = -Z*Z/real(n*n) -0.1;
    // orb_E[j] = 1.5 + 2*(n-1) - l + 0.1;
    orb_w[j] = MAX(0, MIN(2*(2*l+1),Z-nfill));

    nfill += orb_w[j]; j++;
  }

  // show initial electronic configuration
  printf("\nInitial e- config of Z=%i\n",Z);
  for(i = 0; i < nlvl; i++)
    printf("n=%i l=%i  #=%i  %20.15f\n", orb_n[i], orb_l[i], orb_w[i], orb_E[i]);
  printf("\n");

  for(int scfit = 0; scfit < 1; scfit++) {

    // clear total charge density
    for(i = 0; i < N; i++) si[i] = 0.0;

    // find exact energy levels and radial wave functions
    // solving for [-(d/dr)^2 + l(l+1)/r^2 + V(r)]P_nl(r) = E_nl P_nl(r)
    // where P_nl(r) = r R_nl(r) 
    for(i = 0; i < nlvl; i++) {
  
      int it = 0;
      real Emin = -1.e20, Emax = 0.0;
  //    real Emin = 0.0, Emax = 100.0;
      real& E = orb_E[i];
      n = orb_n[i]; l = orb_l[i];
      for(j = 0; j < N; j++) orb_S[i][j] = 0.0;
  
      // if( orb_n[i] != 3 || orb_l[i] != 2 ) continue;
  
      beginning:
  
      // E = it*0.01;
      if( ++it > 1000 )
        { printf("# exceeded max iterations\n"); continue; return -1; }
  
      // check if E is bounded
      bool boundsonE = false;
      if( Emin > -0.99e20 && Emax < 0.0 ) 
  //    if( Emin > -1.0 && Emax < 0.99e20 ) 
        boundsonE = true;
  
      // find closest grid point to classical turning point
      int TP; real C = l*(l+1.0);
      for(TP = N-1; TP > 0 && E < (rV[TP] + C/r[TP])/r[TP]; TP--);  
        if(TP == 0) { Emin = E; E = 0.75*E; goto beginning; }
  //    if(TP == 0) { Emin = E; E = 1.25*E; goto beginning; }
      // TP = 8000;
  
      // printf("# n=%i l=%i  E=%20.15f  TP=%i  it=%i\n",
      //  orb_n[i], orb_l[i], orb_E[i], TP, it);
   
      // perform 4-step Adams-Bashforth-Moulton integration
      // change independent variable to x = log(r/mu):
      //   [-(d/dx)^2 + d/dx + l(l+1) + r^2 (V-E)] P = 0
      // define Q = dP/dx to change to first order equation:
      //   d[P, Q]/dx = [Q, fP+Q] 
      // where f(x) = l(l+1) + r(x)^2 ( V(x)-E )
  
      const real a[4] =  {55./24.*dx,-59./24.*dx,37./24.*dx,-9./24.*dx};
      const real b[4] =  {3./8.*dx,19./24.*dx,-5./24.*dx,1./24.*dx};
      const int ncorr = 3;
      real fP[4], al[4], c[5], dP, dQ;
      real k1P,k1Q,k2P,k2Q,k3P,k3Q,k4P,k4Q,f23;
  
      // reset no. of nodes
      int nodes = 0, reqnodes = n - l - 1;
      real sign = 1.0;
   
      // precompute f(x), g(x)=f(x)/r(x)^2 for all x
      for(j = 0; j < N; j++) f[j] = C + r[j]*(rV[j]-r[j]*E);
      for(j = 0; j < N; j++) g[j] = f[j]/(r[j]*r[j]);
  
      // define taylor series expansion to O(r^5) about r = 0: 
      //   P(x) = sum(n=0..4, c_n r^n)
      // note: does not use the potential
      c[0] = 0.0; c[1] = 1.0; c[2] = 0.5;
    	c[3] = (1.0+C)/6.0; c[4] = (1.0+2.0*C)/24.0;
  
      enum solvede {TE, RK4, ABM4};
      solvede starteq = RK4;
      solvede solveeq = ABM4;
  
      switch( starteq ) {
  
        case TE: // TODO: Not working?
        // set boundary condition P_nl(r=0) = 0, dP/dr(r=0) = a
        // recall index=0 corresponds to r=r_min, and Q = dP/dx
    		// fill in P and Q for first four indices using taylor
        
        for(j = 0; j < 4; j++) {
    		  P[j] = 0.0; Q[j] = 0.0;
          for(k = 4; k >= 0; k--) {
    			  P[j] =   c[k] + r[j]*P[j];
    				Q[j] = k*c[k] + r[j]*Q[j];
    		  }
          fP[3-j] = f[j]*P[j];
        } break;
  
        case RK4:
        // set boundary condition P_nl(r=0) = 0, dP/dr(r=0) = a
        // recall index=0 corresponds to r=r_min, and Q = dP/dx
        // make sure boundary condition matches taylor expansion
  
        P[0] = 0.0; Q[0] = 0.0; 
        for(k = 4; k >= 0; k--) {
    		  P[0] =   c[k] + r[0]*P[0];
    			Q[0] = k*c[k] + r[0]*Q[0];
    		} 
        fP[3] = f[0]*P[0];
   
        // compute first four values of [P,Q] using RK4
        x = xi;
        for(j = 1; j < 4; j++) {
  
          // get half-point for f (do better?)
          f23 = (f[j-1]+f[j])/2.0;
  
          // advance point using Runge-Kutta scheme
          k1P = dx*Q[j-1];           k1Q = dx*f[j-1]*P[j-1] + k1P;
          k2P = dx*(Q[j-1]+k1Q/2.0); k2Q = dx*f23*(P[j-1]+k1P/2.0) + k2P;
          k3P = dx*(Q[j-1]+k2Q/2.0); k3Q = dx*f23*(P[j-1]+k2P/2.0) + k3P;
          k4P = dx*(Q[j-1]+k3Q);     k4Q = dx*f[j]*(P[j-1]+k3P) + k4P;
          P[j] = P[j-1] + (k1P + 2.0*k2P + 2.0*k3P + k4P)/6.0;
          Q[j] = Q[j-1] + (k1Q + 2.0*k2Q + 2.0*k3Q + k4Q)/6.0;
          fP[3-j] = f[j]*P[j]; x + dx;
  
        } break; 
  
      }
        
      switch( solveeq ) {
  
        case RK4:
  
        x = xi + 3.0*dx;
        for(j = 4; j <= TP; j++) {
  
          // get cubic fit of r V(r) for next four grid points
          /* al[0] = x; al[1] = x+dx; al[2] = x+2.*dx; al[3] = x+3.*dx;
          getcubicfit(al, &f[j-1], c); 
          for(f23 = 0.0, k = 3; k >= 0; k--) 
            f23 = f23*(x+dx/2.0) + c[k]; */ 
  
          // get half-point for f (does better?)
          f23 = (f[j-1]+f[j])/2.0;
  
          // advance point using Runge-Kutta scheme
          k1P = dx*Q[j-1];           k1Q = dx*f[j-1]*P[j-1] + k1P;
          k2P = dx*(Q[j-1]+k1Q/2.0); k2Q = dx*f23*(P[j-1]+k1P/2.0) + k2P;
          k3P = dx*(Q[j-1]+k2Q/2.0); k3Q = dx*f23*(P[j-1]+k2P/2.0) + k3P;
          k4P = dx*(Q[j-1]+k3Q);     k4Q = dx*f[j]*(P[j-1]+k3P) + k4P;
          P[j] = P[j-1] + (k1P + 2.0*k2P + 2.0*k3P + k4P)/6.0;
          Q[j] = Q[j-1] + (k1Q + 2.0*k2Q + 2.0*k3Q + k4Q)/6.0;
          x + dx;
  
          // check for zero crossing
          if( P[j]*sign < 0.0 )
            { nodes++; sign = -sign;}
  
        } break;
  
        case ABM4:
  
        // proceed with 4-step Adams-Bashforth-Moulton integrator
        for(j = 4; j <= TP; j++) {
    
          // predictor
          dP = a[0]*Q[j-1] + a[1]*Q[j-2] + a[2]*Q[j-3] + a[3]*Q[j-4];
          dQ = a[0]*fP[0]  + a[1]*fP[1]  + a[2]*fP[2]  + a[3]*fP[3] + dP;
          P[j] = P[j-1] + dP; Q[j] = Q[j-1] + dQ;
          
          // corrector
          for(k = 0; k < ncorr; k++) {
            dP = b[0]*Q[j]      + b[1]*Q[j-1] + b[2]*Q[j-2] + b[3]*Q[j-3];
            dQ = b[0]*f[j]*P[j] + b[1]*fP[0]  + b[2]*fP[1]  + b[3]*fP[2] + dP;
            P[j] = P[j-1] + dP; Q[j] = Q[j-1] + dQ; 
          }
    
          // update array f P to coincide with last four values
          fP[3] = fP[2]; fP[2] = fP[1]; fP[1] = fP[0]; fP[0] = f[j]*P[j];
      
          // check for zero crossing
          if( P[j]*sign < 0.0 )
            { nodes++; sign = -sign;}
  
        } break;
      }
  
      // check for too many nodes
      if( nodes > reqnodes )
        if( boundsonE == true ) 
          { Emax = E; E = (Emax + Emin)/2.0; goto beginning; }
        else
          { Emax = E; E = 1.25*E; goto beginning; } 
  //        { Emax = E; E = 0.75*E; goto beginning; } 
  
      // check for too few nodes
      if( nodes < reqnodes )
        if( boundsonE == true )
          { Emin = E; E = (Emax + Emin)/2.0; goto beginning; }
        else
          { Emin = E; E = 0.75*E; goto beginning; } 
  //        { Emin = E; E = 1.25*E; goto beginning; } 
      
      // determine suitable point to begin WKB tail
      // borrow rule of thumb from original program
      real WKBrad = ( (orb_n[i] > 1? 8:5) + orb_l[i] )*r[TP];
  //    real WKBrad = 8.0;
      for(j = TP; j < N-1 && r[j] < WKBrad; j++);
      int WKB = j;
  
      // printf("# n=%i l=%i  E=%20.15f  WKB=%i\n",
      //  orb_n[i], orb_l[i], orb_E[i], WKB);
   
      // fill in tail using WKB approximation
      real A = 1.0e-03*sqrt(sqrt(f[WKB]/(r[WKB]*r[WKB])));
      real kappa = 0.0;
      for(j = WKB; j < N; j++) {
  
         kappa += 0.5*(sqrt(g[j])+sqrt(g[j-1]))*(r[j]-r[j-1]);
         if( j == WKB ) kappa = 0.0;
      
         real ep = der(dx, g, j);
  
         P[j] = A*exp(-kappa)/sqrt(sqrt(g[j]));
         Q[j] = -(sqrt(g[j])*r[j]+ep/(4.0*g[i]))*P[j];
         if( j-WKB <= 3 ) fP[j-WKB] = f[j]*P[j];   
      }
  
      // now perform inward integration
      // proceed with Adams-Bashforth-Moulton backwards 
  
      // proceed with 4-step Adams-Bashforth-Moulton integrator
      for(j = WKB-1; j > TP; j--) {
  
        // predictor
        dP = a[0]*Q[j+1] + a[1]*Q[j+2] + a[2]*Q[j+3] + a[3]*Q[j+4];
        dQ = a[0]*fP[0]  + a[1]*fP[1]  + a[2]*fP[2]  + a[3]*fP[3] + dP;
        P[j] = P[j+1] - dP; Q[j] = Q[j+1] - dQ;
        
        // corrector
        for(k = 0; k < ncorr; k++) {
          dP = b[0]*Q[j]      + b[1]*Q[j+1] + b[2]*Q[j+2] + b[3]*Q[j+3];
          dQ = b[0]*f[j]*P[j] + b[1]*fP[0]  + b[2]*fP[1]  + b[3]*fP[2] + dP;
          P[j] = P[j+1] - dP; Q[j] = Q[j+1] - dQ; 
        }
  
        // update array f P to coincide with last four values
        fP[3] = fP[2]; fP[2] = fP[1]; fP[1] = fP[0]; fP[0] = f[j]*P[j];
      }

      // proceed one extra step
      real Pextra, Qextra;
         // predictor
        dP = a[0]*Q[j+1] + a[1]*Q[j+2] + a[2]*Q[j+3] + a[3]*Q[j+4];
        dQ = a[0]*fP[0]  + a[1]*fP[1]  + a[2]*fP[2]  + a[3]*fP[3] + dP;
        Pextra = P[j+1] - dP; Qextra = Q[j+1] - dQ;
        
        // corrector
        for(k = 0; k < ncorr; k++) {
          dP = b[0]*Qextra      + b[1]*Q[j+1] + b[2]*Q[j+2] + b[3]*Q[j+3];
          dQ = b[0]*f[j]*Pextra + b[1]*fP[0]  + b[2]*fP[1]  + b[3]*fP[2] + dP;
          Pextra = P[j+1] - dP; Qextra = Q[j+1] - dQ; 
        }
  
      // calculate logarithmic derivatives on both sides of TP
      real ld1 = lder(dx, P, TP);
      real ld2 = lder2(dx, P, TP);
      // real lld = (ld1 + dx*ld2)/((P[TP] + dx*ld1)*r[TP+1]);
      real lld = ld1/(P[TP]*r[TP]);
      real Psave = P[TP];
      P[TP] = Pextra;
      // real rld = rder(dx, P, TP+1)/(P[TP+1]*r[TP+1]);
      real rld = rder(dx, P, TP)/(P[TP]*r[TP]);
      P[TP] = Psave;
  
      // ensure matching derivatives on both sides
      // todo: use Newton-Raphson or Broyden to improve trial energy
      if( it == 999 ) goto checksolution;

      if( lld-rld > 1.0e-06 )
        if( boundsonE == true )
          { Emin = E; E = (Emax + Emin)/2.0; goto beginning; }
        else
          { Emin = E; E = 0.75*E; goto beginning; } 
  //        { Emin = E; E = 1.25*E; goto beginning; } 
  
      if( lld-rld < -1.0e-06 )
        if( boundsonE == true )
          { Emax = E; E = (Emax + Emin)/2.0; goto beginning; }
        else
          { Emax = E; E = 1.25*E; goto beginning; }
  //        { Emax = E; E = 0.75*E; goto beginning; }
  
      // normalize tail to match outward integration
      checksolution:
      real a0 = (P[TP] + dx*lld*r[TP]*P[TP])/P[TP+1];
      for(j = TP+1; j < N; j++)
        { P[j] = a0*P[j]; Q[j] = a0*Q[j]; }

      // normalize wave function
      a0 = 0.0;
      for(j = 1; j < N; j++)
        a0 += 4.*pi*P[j]*P[j]*(r[j]-r[j-1]);
      for(j = 0; j < N; j++)
        P[j] /= sqrt(a0);

      // double-check solution
      /* real error = 0.;
      for(j = 0; j < N; j++) {
  
        // these formulas have too much error to cross-check 
        // predictor-corrector method
        
//      x = xi + j*dx; 
//      al[0] = x; al[1] = x+dx; al[2] = x+2.*dx; al[3] = x+3.*dx;
//      getcubicfit(al,&P[j],c);
//
//      real d0Pc = c[0] + c[1]*x + c[2]*x*x + c[3]*x*x*x;
//      real d1Pc = c[1] + 2.*c[2]*x + 3.*c[3]*x*x;
//      real d2Pc = 2.*c[2] + 6.*c[3]*x; 
//
//      real d1Pf = (P[j+1]-P[j-1])/(2.*dx);
//      real d2Pf = (P[j+1]-2.*P[j]+P[j-1])/(dx*dx); 
//      
  
        // use sixth order finite difference formula
        // nevertheless we are precision limited
        real d1P, d2P;
        if(3 <= j && j <= N-3) {
            d1P = der(dx, P, j);
            d2P = der2(dx, P, j);
        } else if( j < 3 ) {
            d1P = rder(dx, P, j);
            d2P = rder2(dx, P, j);
        } else if( N-3 < j ) {
            d1P = lder(dx, P, j);
            d2P = lder2(dx, P, j);
        }
     
        
          if(j < 5 || N-j < 5 || abs(j-TP) < 10 || 
            abs(WKB-j) < 10 || j % (N/100) == 0)
          if( j <= TP )
            printf("%6i %20.15f %20.15f %20.15f %20.15f %i\n", j, r[j], P[j],
              (-d2P + d1P + f[j]*P[j])/P[j], lder(dx, P, j)/P[j], it);
          else
            printf("%6i %20.15f %20.15f %20.15f %20.15f %i\n", j, r[j], P[j],
              (-d2P + d1P + f[j]*P[j]), rder(dx, P, j)/P[j], it); 
         
        error += log10(abs((-d2P + d1P + f[j]*P[j])/P[j]));
      } 
      static int cnt = 0;
      printf("# avg error = %10.5f %d\n", error/N, cnt++);
      printf("\n"); */
      // goto beginning;
  
      // calculate 4 pi r^2 rho(r) (i.e. the radial charge density)
      printf("# converged n=%i l=%i  E=%20.15f  TP=%i WKB=%i\n",
        orb_n[i], orb_l[i], orb_E[i], TP, WKB);
 
      for(j = 0; j < N; j++) 
        { orb_S[i][j] = 4.*pi*P[j]*P[j]; si[j] += orb_w[i]*orb_S[i][j]; }
    
    }
  
    for(i = 0; i < nlvl; i++)
      printf("# n = %i, l = %i, E = %20.15f\n", orb_n[i], orb_l[i], orb_E[i]); 
 
    // calculate new radius * potential
    // take Slater exchange Vxc = -3 e^2 ((-3/8pi) * rho)^1/3
    // printf("End of eigenvalue solve: iteration = %d\n",scfit);
    fflush(stdout);

    // This is NOT a self-energy correction.
    // The unenclosed charge still contributes to potential (dQ/R)
    real Qselftot = 0.0;
    for(j = 0; j < N; j++)
      if( j > 0 )
        Qselftot += 0.5*(si[j]/r[j]+si[j-1]/r[j-1])*(r[j]-r[j-1]);
//      else
//        Qselftot += si[j];

    real Qenc = 0.0, Qself = 0.0;
    for(j = 0; j < N; j++) {
      if( j > 0 ) {
        Qenc += 0.5*(si[j]+si[j-1])*(r[j]-r[j-1]);
        Qself += 0.5*(si[j]/r[j]+si[j-1]/r[j-1])*(r[j]-r[j-1]);
      }
      else {
        Qenc += si[j]*r[j];
//        Qself += si[j];
      }

      real alpha = 0.70;
      real rVxc = -6.*alpha*pow(3.0*r[j]*si[j]/(32.*pi*pi),1./3.);
      // printf("%20.15f %20.15f %20.15f %20.15f %20.15f %20.15f\n",r[j],Qenc,rV[j],rVxc,Qself,Qselftot);
      real mix = 0.05;
      rV[j] = mix*(2.*(Qenc+r[j]*(Qselftot-Qself)-Z) + rVxc) + (1.-mix)*rV[j];
      // rV[j] = mix*(2.*(Qenc-Z) + rVxc) + (1.-mix)*rV[j];

      // computer Latter tail correction
      rV[j] = MIN(-2., rV[j]);
    }
    printf("\n"); 
  
  }

  // take out the trash
  delete[] r, rV, si, P, Q, f;
  delete[] orb_n, orb_l, orb_E, orb_w;
  for(i = 0; i < nlvl; i++)
    delete[] orb_S[i], orb_P[i];
  delete[] orb_S, orb_P;

  return 0;
}

