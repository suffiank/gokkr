#include "atom.h"
#include "eqsolver.h"
#include "../../util/src/extractkvp.h"
#include <cstdio>
#include <algorithm>
using namespace std;

dble sign(dble val) {
  return (0.0 < val) - (val < 0.0);
}

void atom_t::configure(std::map<std::string,std::string>& kvp) {

  extractkvp(kvp, "Z", Z, 1);

  int ngrid;
  double rmin, rmax;
  extractkvp(kvp, "N", ngrid, 10000);
  extractkvp(kvp, "tf_rmin", rmin, 0.0001);
  extractkvp(kvp, "tf_rmax", rmax, 200.0);
  set_logarithmic_grid(ngrid, rmin, rmax);

  string key = "verbose";
  verbose = false;
  if( kvp.find(key) != kvp.end() )
    if( kvp[key] == "true" ) 
      verbose = true;
}

void atom_t::set_logarithmic_grid(int N, dble rmin, dble rmax) {

  // setup a logarithmic grid
  // use Thomas-Fermi universal scaling factor: 
  //   r = mu r_TF = (1/2)(3pi/4)^(2/3) Z^(-1/3) r_TF

  this->N = N;
  xi = log(rmin);
  dx = ( log(rmax)-xi )/(N-1);
  xf = xi + (N-1)*dx;  
  mu = 0.88534138*pow(Z,-1.0/3.0);

  r.resize(N);
  rV.resize(N);
  sigma.resize(N);

  dble x = xi;
  for(int i = 0; i < N; i++)
    r[i] = mu*exp(xi + i*dx); 
}

void atom_t::fill_madelung_orbitals() {

  // use Madelung's rule to guess highest filled nl subshell
  nlvl = 0;
  int n, l, nfill = 0;
  for(int i = 1; i <= 12; i++)
  for(n = i/2+1, l = (i-1)/2; n <= i; n++, l--)
  if( nfill < Z ) 
    { nfill += MIN(2*(2*l+1), Z-nfill); nlvl++; }

  // fill in a guess electronic configuration
  // keep one extra level in case filling changes
  // take initial energy guess as from V(r) = -2Z/r`
  nlvl++;
  orbital.resize(nlvl);
  for(int i = 0; i < nlvl; i++) {
    orbital[i].wavef.resize(N);
    orbital[i].sigma.resize(N);
  }

  int j = 0; nfill = 0;
  for(int i = 1; i <= 12; i++)
  for(n = i/2+1, l = (i-1)/2; n <= i; n++, l--)
  if( j < nlvl ) {
    orbital[j].n = n; orbital[j].l = l;
    orbital[j].E = -Z*Z/dble(n*n)-0.01;
    orbital[j].occ = MAX(0, MIN(2*(2*l+1),Z-nfill));

    nfill += orbital[j].occ; j++;
  }
}

void atom_t::fill_thomas_fermi_potential() {

  // fill in Latter's fit of Thomas-Fermi potential
  // see Phys.Rev. v99 p510 by Latter
  // see DFT by Parr & Wang, Sec 6.2 pg 110 (in Hartree units)

  dble rx;
  for(int i = 0; i < N; i++) {
    dble x = r[i]/mu; rx = sqrt(x); 
    rV[i] = 1.0 + 0.02747*rx + 1.243*x - 0.1486*x*rx + 0.2302*x*x + \
            0.007298*x*x*rx + 0.006944*x*x*x;
    rV[i] = 2.0*Z/rV[i];
    sigma[i] = 4.0/(3.0*pi)*rV[i]*sqrt(rV[i]*r[i]);
    rV[i] = -rV[i];
  }
}

void atom_t::fill_coloumb_potential() {
  for(int i = 0; i < N; i++)
    rV[i] = -2.0*Z; 
}

void atom_t::fill_square_well(dble U0, dble R0) {
  for(int i = 0; i < N; i++)
    rV[i] = (r[i] < R0? r[i]*U0 : 0.0); 
}

void atom_t::fill_harmonic_oscillator(dble omega) {
  for(int i = 0; i < N; i++)
    rV[i] = 0.25*omega*omega*r[i]*r[i]*r[i]; 
}
 
int atom_t::recompute_orbital(int lvl) {

  // find exact energy level and radial wave function
  //   for [-(d/dr)^2 + l(l+1)/r^2 + V(r)]P_nl(r) = E_nl P_nl(r)
  // where P_nl(r) = r R_nl(r)
  
  // on logarithmic grid,
  // change independent variable to x = log(r/mu):
  //   [-(d/dx)^2 + d/dx + l(l+1) + r^2 (V-E)] P = 0

  // define Q = dP/dx to change to first order equation:
  //   d[P, Q]/dx = [Q, fP+Q] 
  // where f(x) = l(l+1) + r(x)^2 ( V(x)-E )
 
  int n = orbital[lvl].n;
  int l = orbital[lvl].l;
  int occ = orbital[lvl].occ;
  dble& E = orbital[lvl].E;
  xfunc& P = orbital[lvl].wavef;
  xfunc Q(N);

  int it = 0, i, j, k;
  dble Emin = -1.e20, Emax = 1.e20;
  dble raise = 1.0 + 0.25*sign(E);
  dble lower = 1.0 - 0.25*sign(E);
  dble Eprev = E, Fprev = 1.0;
  bool boundsonE = false;

  // clear radial charge density
  for(i = 0; i < N; i++) 
    orbital[lvl].sigma[i] = 0.0;

  const int maxit = 500;

  beginning:
  
  if( ++it > maxit ) return -it;
  
  // check if E is bounded
  if( Emin > -0.99e20 && Emax < 0.99e20 ) 
    boundsonE = true;
  
  // find closest grid point to classical turning point
  int TP; dble C = l*(l+1.0);
  for(TP = N-1; TP > 0 && E < (rV[TP] + C/r[TP])/r[TP]; TP--);  
  if(TP == 0)   { Emin = E; E = raise*E; goto beginning; }
  if(TP == N-1) { Emax = E; E = lower*E; goto beginning; }

  // select a match point
  int MP = MIN(TP,N-4);

  // reset no. of nodes
  int nodes = 0, reqnodes = n - l - 1;
  dble sign = 1.0;
   
  // precompute f(x), g(x)=f(x)/r(x)^2 for all x
  xfunc f(N), g(N);
  for(i = 0; i < N; i++) { 
    f[i] = C + r[i]*(rV[i]-r[i]*E);
    g[i] = f[i]/(r[i]*r[i]);
  }
  
  // set solver for schrodinger equation 
  enum solvede {TE, RK4, ABM4};
  solvede starteq = TE;
  solvede solveeq = ABM4;

  // find terms in series expansion
  //   f(r) = f0 + f1 r + f2 r^2 + ...
  // cubic fit precise only for r << 1
  dble a[4]; getcubicfit(r, rV, a);
  dble f0 = C, f1 = a[0], f2 = -E+a[1];
 
  // set boundary condition P_nl(r << 1) = A r^(l+1)
  // recall index=0 corresponds to r=r_min, and Q = dP/dx
  interval_t origin(0,3);
  solve_sch_taylor_origin(r, P, Q, l, f1, f2, origin); 
  switch( starteq ) {
    case TE: break;
    case RK4:
      origin.a = 1;
      solve_sch_rk4(dx, P, Q, f, origin);
      break;
  }
  
  // solve schrodinger equation by integrating outward
  interval_t inside(4, MP);
  switch( solveeq ) { 
    case RK4:  nodes = solve_sch_rk4(dx, P, Q, f, inside); break;
    case ABM4: nodes = solve_sch_abm4(dx, P, Q, f, inside); break;
  }
  
  // check for too many nodes
  if( nodes > reqnodes )
    if( boundsonE == true ) 
      { Emax = E; E = (Emax + Emin)/2.0; goto beginning; }
    else
      { Emax = E; E = lower*E; goto beginning; } 
  
  // check for too few nodes
  if( nodes < reqnodes )
    if( boundsonE == true )
      { Emin = E; E = (Emax + Emin)/2.0; goto beginning; }
    else
      { Emin = E; E = raise*E; goto beginning; } 
  
  // calculate logarithmic derivative from the inside
  dble inP = P[MP];
  dble lld = Q[MP]/(P[MP]*r[MP]);
  
  // determine suitable point to begin WKB tail
  // note: borrowed rule of thumb from original program
  dble WKBrad = ( (orbital[lvl].n > 1? 8:5) + orbital[lvl].l )*r[TP];
  for(i = TP; i < N-4 && r[i] < WKBrad; i++);
  int WKB = MIN(MAX(i,MP+4),N-4);
  
  // printf("# n=%i l=%i E=%20.15f TP=%i WKB=%i\n", n, l, E, TP, WKB);
   
  // fill in tail using WKB approximation
  interval_t tail(WKB,N-1);
  solve_sch_wkb_tail(dx, r, P, Q, g, tail);

  // now perform inward integration
  interval_t outside(WKB-1,MP);
  switch( solveeq ) { 
    case RK4:  solve_sch_rk4<-1>(dx, P, Q, f, outside); break;
    case ABM4: solve_sch_abm4<-1>(dx, P, Q, f, outside); break;
  }
  
  // calculate logarithmic derivatives on both sides of MP
  dble rld = Q[MP]/(P[MP]*r[MP]);
  
  // ensure matching derivatives on both sides
  if( lld-rld > 1.0e-10 )
    if( boundsonE == true )
      { Emin = E; E = (Emax + Emin)/2.0; goto beginning; }
    else
      { Emin = E; E = raise*E; goto beginning; } 
  
  if( lld-rld < -1.0e-10 )
    if( boundsonE == true )
      { Emax = E; E = (Emax + Emin)/2.0; goto beginning; }
    else
      { Emax = E; E = lower*E; goto beginning; }

  // new: use secant method to guess new trial energy
  /* dble dFdE = (lld-rld-Fprev)/(E-Eprev);
  Eprev = E; Fprev = lld-rld;
  if( abs(lld-rld) > 1.0e-10 && boundsonE == true ) {
    dble Esecant = E - (lld-rld)/dFdE;
    if( Emin < Esecant && Esecant < Emax )
      E = Esecant;
    else
      E = (Emax + Emin)/2.0;
    goto beginning;
  } */
  
  // ensure continuity of solutions at matching point
  dble a0 = inP/P[MP];
  for(i = MP-3; i < N; i++)
    { P[i] = a0*P[i]; Q[i] = a0*Q[i]; }

  // normalize wave function: integral( 4 pi r^2 R^2 dr ) = 1
  // note: need to perform more accurate integration
  a0 = 0.0;
  for(i = 1; i < N; i++)
    a0 += 4.*pi*P[i]*P[i]*(r[i]-r[i-1]);
  for(i = 0; i < N; i++)
    P[i] /= sqrt(a0);

  // fill in radial charge density
  for(i = 0; i < N; i++) 
    orbital[lvl].sigma[i] = 4.*pi*P[i]*P[i]; 

  // return number of iterations required
  if( verbose )
    printf("# n = %i, l = %i, occ = %2i, E = %22.15f, it = %i\n", n, l, occ, E, it);
  return it;
}

void atom_t::recompute_charge() {

  // clear total radial charge density
  for(int i = 0; i < N; i++) 
    sigma[i] = 0.0;

  // add each orbital radial charge density to total
  for(int i = 0; i < nlvl; i++)
    for(int j = 0; j < N; j++) 
      sigma[j] += orbital[i].occ * orbital[i].sigma[j];
}

void atom_t::recompute_potential() {

  // integrate to find charge enclosed by radius r as well as
  // also, potential from charges outside r is Vout(dQ(R)) = dQ(R)/R
  // warning: integration may be of low order
  xfunc Qenc(N), Vout(N); 
  Qenc[0] = sigma[0]*r[0]; 
  Vout[0] = sigma[0];
  for(int i = 1; i < N; i++) {
    Qenc[i] = Qenc[i-1] + 0.5*(sigma[i]+sigma[i-1])*(r[i]-r[i-1]);
    Vout[i] = Vout[i-1] + 0.5*(sigma[i]/r[i]+sigma[i-1]/r[i-1])*(r[i]-r[i-1]);
  }

  // compute Slater exchange potential
  dble alpha = 0.70;
  xfunc rVxc(N);
  for(int i = 0; i < N; i++)
    rVxc[i] = -6.*alpha*pow(3.0*r[i]*sigma[i]/(32.*pi*pi),1./3.);
 
  // combine potentials into single mean-field potential
  // include Latter tail correction (i.e. force fall off to -2/r)
  for(int i = 0; i < N; i++) {
    rV[i] = 2.*(Qenc[i]-Z) + 2.*r[i]*(Vout[N-1]-Vout[i]) + rVxc[i];
    rV[i] = MIN(-2., rV[i]);
  }
 
}

void atom_t::solve_orbitals() {
  // solve for energy levels and radial wave functions
  for(int i = 0; i < nlvl; i++)
    recompute_orbital(i);
}

void atom_t::sort_orbitals() {
  for(int i = 0; i < nlvl; i++)
    for(int j = i; j > 0 && orbital[j].E < orbital[j-1].E; j--)
      swap(orbital[j-1], orbital[j]);
}

int atom_t::solve_self_consistent() {

  const int maxit = 100;

  for(int it = 0; it < maxit; it++) {
    printf("# scf iteration = %d\n",it);

    // solve for energy levels and radial wave functions
    #pragma omp parallel for
    for(int i = 0; i < nlvl; i++)
      recompute_orbital(i);
    
    // calculate total radial charge density
    recompute_charge();

    // solve for a new potential
    xfunc rVcurr = rV;
    recompute_potential();

    // only retain a fraction of new potential
    dble mix = 0.10;
    for(int i = 0; i < N; i++)
      rV[i] = (1.0-mix)*rVcurr[i] + mix*rV[i];
  }

  return 0;
}

void atom_t::solve_reg_fn(int l, dble rad, cplx E, zfunc *Rl, cplx *tl, dble *dl) {

  zfunc P(N), Q(N), f(N);

  // precompute f(x) = l(l+1) + r(x)^2 ( V(x)-E )
  cplx C = l*(l+1.0);
  for(int i = 0; i < N; i++)  
    f[i] = C + r[i]*(rV[i]-r[i]*E);

  // find approximate position where radius cuts off
  // warning: this will lead to numerical inaccuracy
  int MP = 0;
  for(int i = 0; i < N && r[i] < rad; i++) MP = i;

  // warning: for debugging purposes only!
  // this is meant for finite square well
  // MP = 8121;
  // f[MP] = C + r[MP]*r[MP]*(-10.0-E);

  // find terms in series expansion
  //   f(r) = f0 + f1 r + f2 r^2 + ...
  // cubic fit precise only for r << 1
  cplx a[4]; getcubicfit(r, rV, a);
  cplx f0 = C, f1 = a[0], f2 = -E+a[1];
 
  // set boundary condition P_nl(r << 1) = A r^(l+1)
  // recall index=0 corresponds to r=r_min, and Q = dP/dx
  interval_t origin(0,3);
  solve_sch_riccati<bessel_j>(E, l, r, P, Q, origin); 
 
  // integrate outward from spherical bessel at origin
  interval_t nonzero(4, MP);
  solve_sch_abm4(dx, P, Q, f, nonzero);

  // matching to outside free solution to get t-matrix
  // see Zabloudil (pg. 63)
  cplx p = sqrt(E);
  if( p.imag() < 0.0 ) p = -p;

  cplx ld  = (Q[MP]-P[MP])/(P[MP]*r[MP]);
  cplx jl  = sph_bessel_jl(l, p*r[MP]);
  cplx nl  = sph_bessel_yl(l, p*r[MP]);
  cplx jl1 = sph_bessel_jl(l+1, p*r[MP]);
  cplx nl1 = sph_bessel_yl(l+1, p*r[MP]);

  cplx jlp = -jl1 + double(l)/(p*r[MP])*jl;
  cplx nlp = -nl1 + double(l)/(p*r[MP])*nl;

  cplx cdl = atan( (ld*jl-p*jlp)/(ld*nl-p*nlp) );
  *tl = -sin(cdl)*exp(im*cdl)/p;

  // change jost function into regular solution:
  //   Rl = jl - ip tl hl for r > R_cs
  // also, phase dl = phase( limit( r->0, Rl(r)/jl(r) ) )
  
  cplx hl = jl + im*nl;
  cplx A  = r[MP]*( jl - im*p*(*tl)*hl );
 
  *dl = imag_part( log(A/P[MP]) );

  for(int i = 0; i <= MP; i++)
    { P[i] = A/P[MP] * P[i]; Q[i] = A/P[MP] * Q[i]; }

  for(int i = 0; i <= MP; i++)
    (*Rl)[i] = P[i]/r[i];
  for(int i = MP; i < N; i++) {
    jl = sph_bessel_jl(l, p*r[i]);
    hl = sph_bessel_hl1(l, p*r[i]);
    (*Rl)[i] = jl - im*p*hl*(*tl);
  }
}

void atom_t::solve_irr_fn(int l, dble rad, cplx E, zfunc *Sl) {

  zfunc P(N), Q(N), f(N);

  // precompute f(x) = l(l+1) + r(x)^2 ( V(x)-E )
  cplx C = l*(l+1.0);
  for(int i = 0; i < N; i++)  
    f[i] = C + r[i]*(rV[i]-r[i]*E);

  // find approximate position where radius cuts off
  // warning: this will lead to numerical inaccuracy
  int MP = 0;
  for(int i = 0; i < N && r[i] < rad; i++) MP = i;

  // fill in -ip hl+(pr) outside radius
  interval_t zeropot(MP, N);
  solve_sch_riccati<bessel_h1>(E, l, r, P, Q, zeropot);

  cplx p = sqrt(E);
  if( p.imag() < 0.0 ) p = -p;
  for(int i = MP; i < N; i++)
    { P[i] = -im*p*P[i]; Q[i] = -im*p*Q[i]; }

  // integrate inward from outside hankel function
  interval_t nonzero(MP,0);
  solve_sch_abm4<-1>(dx, P, Q, f, nonzero);

  // fill in irregular solution
  for(int i = 0; i < N; i++)
    (*Sl)[i] = P[i]/r[i];
}

/* cplx atom_t::greenfn(cplx E, const vec3& x1, const vec3& x2) {

  zfunc Rl(N), Sl(N);

  cplx G = 0.0;
  for(int l = 0; l <= maxl; l++) {

    // warning: recompute each time?? better to cache?
    cplx tl; dble dl;
    solve_reg_fn(l, rad, E, &Rl, &tl, &dl);
    solve_irr_fn(l, rad, E, &Sl);

    dble rmin = min(x1.mag(), x2.mag());
    dble rmax = max(x1.mag(), x2.mag());
    G += Rl(rmin) * Sl(rmax);
  }

  return G;
} */
