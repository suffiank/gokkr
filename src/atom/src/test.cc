#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

#include "../../util/src/mathfn.h"
#include "atom.h"

char speclett(int l) {
  switch(l) {
    case 0: return 's';
    case 1: return 'p';
    case 2: return 'd';
    case 3: return 'f';
    default: 
      return char(int('g')+l-4);
  };
}

int main(int argc, char *argv[]) {

  atom_t atom;
  int n, l, reqit[10][10];

  // fill an initial set of energy levels
  atom.Z = 100;
  atom.set_logarithmic_grid(10000, 0.0001, 200.0);
  atom.fill_madelung_orbitals();
  goto squarewell;

  // solve for Coloumb potential
  /* atom.Z = 10;
  atom.fill_coloumb_potential();
  for(int i = 0; i < atom.nlvl; i++) {
    n = atom.orbital[i].n; l = atom.orbital[i].l;
    reqit[n][l] = atom.recompute_orbital(i);
  }
  atom.sort_orbitals();

  // print comparison with exact Coloumb energies
  printf("  -- Energy spectra of Columb potential (-20/r) --\n");
  printf("  it   el-     energy-solved        energy-exact\n");
  for(int i = 0; i < atom.nlvl; i++) {
    int n = atom.orbital[i].n, l = atom.orbital[i].l;
    dble Esolve = atom.orbital[i].E;
    dble Eexact = -dble(atom.Z)*atom.Z/(n*n);
    printf("%4d %3d%1c %20.15f %20.15f\n", reqit[n][l], n, speclett(l), 
        Esolve, Eexact); 
  }

  // solve for harmonic oscillator
  atom.fill_harmonic_oscillator(2.0);
  for(int i = 0; i < atom.nlvl; i++) {
    atom.orbital[i].E = +0.1;
    n = atom.orbital[i].n; l = atom.orbital[i].l;
    reqit[n][l] = atom.recompute_orbital(i);
  }
  atom.sort_orbitals();

  // print comparison with exact oscillator energies
  printf("  -- Energy spectra of harm. oscillator  (r*r)  --\n");
  printf("  it   el-     energy-solved        energy-exact\n");
  for(int i = 0; i < atom.nlvl; i++) {
    int n = atom.orbital[i].n, l = atom.orbital[i].l;
    dble Esolve = atom.orbital[i].E;
    dble Eexact = (1.5+2.*(n-1-l)+l)*2.0;
    printf("%4d %3d%1c %20.15f %20.15f\n", reqit[n][l], n, speclett(l), 
        Esolve, Eexact); 
  }
  // exit(0);
  */
 
  // solve for finite square well
  // set r[8121] = 2.5 exactly
  // note: you can then fix the match point to R = 2.5 exactly
  // you also have to be careful to evaluate only left & right deriv.
  // AND to watch out for the corrector part of ABM4 (use RK4)
  squarewell:
  atom.Z = 100;
  atom.set_logarithmic_grid(10000, 0.0001, 199.953337378160);
 
  // calculate orbitals
  atom.Z = 130;
  atom.fill_madelung_orbitals();
  atom.fill_square_well(-10.0, 2.5);
  /* for(int i = 0; i < atom.nlvl; i++) {
    atom.orbital[i].E = -0.5;
    n = atom.orbital[i].n; l = atom.orbital[i].l;
    reqit[n][l] = atom.recompute_orbital(i);
  }
  atom.sort_orbitals(); 
 
  printf("  -- Energy spectra of square well (U=10,R=2.5) --\n");
  printf("   note: discontinuity in V(r) will affect solver\n");
  printf("  it   el-     energy-solved        energy-exact\n");
  for(int i = 0; i < atom.nlvl; i++) {
    n = atom.orbital[i].n, l = atom.orbital[i].l;
    dble Esolve = atom.orbital[i].E;
    dble Eexact = 1.0;
    
    if(n == 1 && l == 0) Eexact = -8.76168059068852;
    if(n == 2 && l == 1) Eexact = -7.47828564237501;
    if(n == 3 && l == 2) Eexact = -5.87431062290218;
    if(n == 2 && l == 0) Eexact = -5.13934477976623; 
    if(n == 4 && l == 3) Eexact = -3.97633438582846; 
    if(n == 3 && l == 1) Eexact = -2.77725285523125; 
    if(n == 5 && l == 4) Eexact = -1.81241558709838; 
    if(n == 4 && l == 2) Eexact = -0.29747036554810; 
    if(n == 3 && l == 0) Eexact = -0.01938732849664; 
    
    if( Eexact < 0.0 ) 
      printf("%4d %3d%1c %20.15f %20.15f\n", reqit[n][l], n, speclett(l), 
          Esolve, Eexact);
    else if( reqit[n][l] > 0 )
        printf("%4d %3d%1c %20.15f %20s\n", 
            reqit[n][l], n, speclett(l), Esolve, "no solution");
  } */

  // printf("energy, phase-shift\n");
  const int maxl = 7;
  dble R = 2.5;
  cplx d[maxl], _d[maxl], _d0e;
  cplx t[maxl], _t[maxl];

  for(int i = 0; i < maxl; i++)
    { _d[i] = _t[i] = 0.0; }
  _d0e = _d[0]*im;
  for(dble E = -10.1; E < 50.0; E += 0.001) { 

    // exact s-wave phase shift
    cplx k = sqrt( cplx(E,0.001) );
    if( k.imag() < 0.0 ) k = -k;
    cplx k0 = sqrt(10.0);
    cplx K = sqrt(k*k + 10.0);
    if( K.imag() < 0.0 ) K = -K;
    // cplx d0exact = -k*R + atan(k*tan(K*R)/K);
    cplx z = k/K * tan(K*R);
    cplx d0exact = -k*R + log( (1.0+im*z)/(1.0-im*z) )/(2.*im); 
    cplx tmat = log(-k*sin(d0exact)*exp(im*d0exact));
    while( imag_part(d0exact-_d0e) > pi )  d0exact -= 2.0*pi*im; 
    while( imag_part(d0exact-_d0e) < -pi ) d0exact += 2.0*pi*im;
    _d0e = d0exact;
    printf("%20.15f %20.15f",E,d0exact.real());
 
    // numerical phase shift
    for(int i = 0; i < maxl; i++) {
      dble dl; cplx tl; zfunc Rl(10000);
      atom.solve_reg_fn(i, 2.5, cplx(E,0.001), &Rl, &tl, &dl);
      d[i] = dl; t[i] = log(tl)/im;
      // d[i] = atom.phase_shift(i, 2.5, cplx(E,0.001));

      while( real_part(d[i]-_d[i]) > pi )  d[i] -= 2.0*pi; 
      while( real_part(d[i]-_d[i]) < -pi ) d[i] += 2.0*pi;
      while( real_part(t[i]-_t[i]) > pi )  t[i] -= 2.0*pi; 
      while( real_part(t[i]-_t[i]) < -pi ) t[i] += 2.0*pi;
      printf(" %20.15f %20.15f",d[i].real(),t[i].real());
      _d[i] = d[i]; _t[i] = t[i];
    }
    printf("\n");
  }
}
