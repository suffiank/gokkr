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

  // solve for finite square well
  // set r[8121] = 2.5 exactly
  // note: you can then fix the match point to R = 2.5 exactly
  // you also have to be careful to evaluate only left & right deriv.
  // AND to watch out for the corrector part of ABM4 (use RK4)
  atom.Z = 100;
  atom.set_logarithmic_grid(10000, 0.0001, 199.953337378160);

 
  atom.Z = 130;
  atom.fill_madelung_orbitals();
  atom.fill_square_well(-10.0, 2.5);
  const int maxl = 7;
  dble R = 2.5;
  cplx d[maxl], _d[maxl], _d0e;
  cplx t[maxl], _t[maxl];
  
  printf("# %18s %20s %20s %20s %20s %20s %20s %20s\n",
      "energy","l=0","l=1","l=2","l=3","l=4","l=5","l=6");

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
    // printf("%20.15f %20.15f",E,d0exact.real());
    printf("%20.15f",E);
 
    // numerical phase shift
    for(int i = 0; i < maxl; i++) {
      dble dl; cplx tl; zfunc Rl(10000);
      atom.solve_reg_fn(i, 2.5, cplx(E,0.001), &Rl, &tl, &dl);
      d[i] = dl; t[i] = log(tl)/im;

      while( real_part(d[i]-_d[i]) > pi )  d[i] -= 2.0*pi; 
      while( real_part(d[i]-_d[i]) < -pi ) d[i] += 2.0*pi;
      while( real_part(t[i]-_t[i]) > pi )  t[i] -= 2.0*pi; 
      while( real_part(t[i]-_t[i]) < -pi ) t[i] += 2.0*pi;
      printf(" %20.15f",d[i].real());
      _d[i] = d[i]; _t[i] = t[i];
    }
    printf("\n");
  }
}
