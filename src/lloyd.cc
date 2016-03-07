#include <iostream>
using namespace std;

#include "util/src/loadkvp.h"
#include "util/src/mat3.h"
#include "util/src/mathfn.h"
#include "util/src/timer.h"
#include "crystal/src/crystal.h"
#include <Eigen/Dense>
using namespace Eigen;

void track_phase(cplx& phi, cplx& prev) {
  while( imag_part(phi-prev) > pi  ) phi -= 2.0*pi*im; 
  while( imag_part(phi-prev) < -pi ) phi += 2.0*pi*im;
  prev = phi;
}

int main(int argc, char **argv) {

  map<string,string> kvp;
  bool good = safeloadkvp("input.txt",&kvp);
  if( !good ) return -1;

  crystal_t crystal(kvp);

  // setup atom grid
  atom_t& atom = crystal.atom[0];
  atom.Z = 100;
  atom.set_logarithmic_grid(10000, 0.0001, 199.953337378160);
 
  // calculate orbitals
  atom.Z = 130;
  atom.fill_madelung_orbitals();
  atom.fill_square_well(-10.0, 2.5);
 
  // make room for tau matrix and previous phases
  int matdim = crystal.green0.get_g0dim();
  cplx *gmat = new cplx[matdim*matdim];
  cplx *pphi = new cplx[crystal.specialk.numk()*matdim];
  cplx *pevl = new cplx[crystal.specialk.numk()*matdim];
  for(int i = 0; i < matdim; i++) pphi[i] = 0.0;
  for(int i = 0; i < matdim; i++) pevl[i] = 0.0;

  // map a matrix from eigen to the above space
  typedef Matrix<cplx, Dynamic, Dynamic> MatrixXz;
  Map< MatrixXz > gmatE(gmat, matdim, matdim);
  Matrix<cplx,Dynamic,1> evals(matdim);

  const int maxl = 3;
  double _dl[maxl+1];
  for(int l = 0; l <= maxl; l++)
    _dl[l] = 0.0;

  printf("nspecialk = %d\n",crystal.specialk.numk());
  for(cplx E = cplx(-10.0,0.0005); E.real() < 5.0; E += 0.0003) {

    printf("%20.15f",E.real());
    // calculate single-site t-matrix
    cplx tl[maxl+1]; double dl[maxl+1];
    for(int l = 0; l <= maxl; l++) {
      zfunc Rl( atom.N );
      atom.solve_reg_fn(l, 2.5, E, &Rl, &tl[l], &dl[l]);
      cplx zdl = cplx(0.0,dl[l]);
      cplx _zdl = cplx(0.0,_dl[l]);
      track_phase(zdl,_zdl);
      dl[l] = zdl.imag(); _dl[l] = _zdl.imag();
      printf("%20.15f",dl[l]);
    }
 
    // for each special k point
    cplx lloyd = 0.0; dble totalk = 0.0;
    vec3 kvec = 0.5*crystal.bvec[0];
    for(int k = 0; k < crystal.specialk.numk(); k++) {
      totalk += (double)crystal.specialk[k].weight();
      //printf("\nk = %20.15f %20.15f %20.15f\n",
      //    crystal.specialk[k].k.x,crystal.specialk[k].k.y,
      //    crystal.specialk[k].k.z);
  
      // fill in g0 matrix
      // crystal.green0.calc_g0mat(gmat, E, crystal.specialk[k].k);
      crystal.green0.calc_g0mat(gmat, E, kvec);
      for(int i = 0; i < matdim; i++)
        gmatE(i,i) += im*sqrt(E);
 
      // multiply by t-matrix
      for(int l = 0, L = 0; l <= maxl; l++)
      for(int m = -l; m <= l; m++, L++)
        for(int L2 = 0; L2 < matdim; L2++)
          gmat[L*matdim+L2] *= tl[l]; 

      // change to t^-1 - g0
      /* gmatE = -gmatE;
      for(int l = 0, L = 0; l <= maxl; l++)
      for(int m = -l; m <= l; m++, L++)
          gmat[L*matdim+L] += 1.0/tl[l]; */
 
      // change to 1 - t g0
      gmatE = MatrixXz::Identity(matdim,matdim) - gmatE;
      evals = gmatE.eigenvalues();

      /* printf("\n");
      for(int i = 0; i < matdim; i++)
        printf("%20.15f %20.15f\n",evals[i].real(),evals[i].imag());
      printf("\n"); */

      // sort eigenvalues to ensure continuity
      if( E.real() > -0.7000001 ) {
        vector<bool> marked(matdim,false);
        vector<int> perm(matdim,-1);
        for(int i = 0; i < matdim; i++) {
          double min = 1.e10; int minj = -1;
          for(int j = 0; j < matdim; j++)
            if( !marked[j] && abs(evals[j]-pevl[matdim*k+i]) < min ) 
              { min = abs(evals[j]-pevl[matdim*k+i]); minj = j; }
          perm[i] = minj; marked[minj] = true;
          pevl[matdim*k+i] = evals[minj];
        }
      }
      else
        for(int i = 0; i < matdim; i++)
          pevl[matdim*k+i] = evals[i];
  
      // place eigenvalue phase on correct branch
      for(int n = 0; n < matdim; n++) {
        cplx phi = log(pevl[n]);
        track_phase( phi, pphi[k*matdim+n] ); 
        printf(" %20.15f %20.15f",pevl[n].real(),pevl[n].imag());
        // printf("%20.15f %20.15f",phi.real(),phi.imag());
        lloyd += phi * (double)crystal.specialk[k].weight();
      } 
    } 
    lloyd /= -pi*totalk;

    printf(" %20.15f\n",lloyd.imag());
  }

  // free matrix and diagonal phases
  delete [] gmat;
  delete [] pphi;
  delete [] pevl;

}
