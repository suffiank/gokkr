/* written by Suffian Khan Aug 2014 */

#include "green0.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <ctime>
#include <map>

#define USE_MATH_DEFINES
#include <cmath>

#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_integration.h"
#include "ylm.h"

using namespace std;

// fill in gaunt coefficients:
//   integral Y_L1(O) Y_L2(O) conj(Y_Lp(O)) dO
// for O = Omega
// warning: modified from above to agree with MECCA
void green0_t::calc_all_gaunt() {
  try{ 
    if(log_times) timer.begin("Generate Gaunt coefficients");
    fflush(stdout);
    calc_all_gaunt_fixed();
    if(log_times) timer.end();
  } 
  catch(bad_alloc& ba) {
    throw("calc_all_gaunt: insufficient memory\n");
  }
}

// use fixed Gauss-Legendre grid to perform integration
void green0_t::calc_all_gaunt_fixed() {

  // clear any existing gaunt table
  gaunt.clear(); gLpval.clear(); gindex.clear();

  // integration error tolerances
  const double epsabs = 1e-12; 
  double I; 

  // setup gauss-legendre abscissae and weights
  // note: must specify power of two to get machine
  //  precision numbers from the gsl library
  int numw = 1 << (int)ceil(log2((double)numLp));
  bool isodd = numw % 2;
  int nstored = isodd? numw/2+1 : numw/2; 
  gsl_integration_glfixed_table *gltable = 
     gsl_integration_glfixed_table_alloc(numw);
 
  // setup associated legendre polynomials
  // note: gsl only stores x >= 0 abscissae

  // for every absiccae
  vector< vector<double> > sphL(numw);
  for(int i = 0; i < numw; i++)
    sphL[i].resize(numLp); 

  for(int i = 0, j = 0; i < nstored; i++) {

    // get associated legendres at x 
    double x = gltable->x[i];
    calc_plm(maxlp, &sphL[j][0],x,sqrt(1.0-x*x));

    // apply spherical harmonic normalization
    // note: l increments before m(>=0)
    for(int k = 0; k < numLp; k++) 
      sphL[j][k] *= Alm[k];
    j++;

    // if x = 0, go to next abscissa
    if(isodd && i==0) continue;
  
    // otherwise include negative abscissa
    calc_plm(maxlp, &sphL[j][0],-x,sqrt(1.0-x*x));
    for(int k = 0; k < numLp; k++) 
      sphL[j][k] *= Alm[k];
    j++;    
  }  
 
  // go through every (L1,L2) pair
  for(int l1 = 0, L1 = 0; l1 <= maxl; l1++)
  for(int m1 = -l1; m1 <= l1; m1++, L1++)
  for(int l2 = 0, L2 = 0; l2 <= maxl; l2++)
  for(int m2 = -l2; m2 <= l2; m2++, L2++) {
  
    // set pointer into gaunt array for (L1,L2)
    gindex.push_back(gaunt.size());

    // lp must satisfy the triangle inequality
    int lp_lo = abs(l1-l2), lp_hi = l1+l2;
    for(int lp = lp_lo; lp <= lp_hi; lp++)
    for(int mp = -lp; mp <= lp; mp++) {
    
      // gaunt is zero if m values don't match
      // if( m1 + m2 != mp ) continue; // <-- correct
      if( m2 != m1 + mp ) continue; // <-- agree w/ MECCA?

      // gaunt is zero if overall parity odd 
      if( (l1+l2+lp) % 2 == 1 ) continue;

      // indexing increments l(>=m) before m(=>0)
      int k1 = abs(m1)*maxlp-abs(m1)*(abs(m1)-1)/2+l1;
      int k2 = abs(m2)*maxlp-abs(m2)*(abs(m2)-1)/2+l2;
      int kp = abs(mp)*maxlp-abs(mp)*(abs(mp)-1)/2+lp;

      // perform gauss-legendre non-adaptive integration
      I = 0.0;
      for(int i = 0, j = 0; i < nstored; i++) {

        double x = gltable->x[i], w = gltable->w[i];
        I += w*sphL[j][k1]*sphL[j][k2]*sphL[j][kp]; j++;
        if(isodd && i == 0) continue;
        I += w*sphL[j][k1]*sphL[j][k2]*sphL[j][kp]; j++;
      }
      I *= tpi;

      // minus factors since gsl harmonics only for m > 0
      double sign = 1.0;
      if( m1 >= 0 && m1%2 == 1 ) sign = -sign;
      if( m2 >= 0 && m2%2 == 1 ) sign = -sign;
      if( mp >= 0 && mp%2 == 1 ) sign = -sign;
      I = sign*I;

      // append nonzero gaunt value to list
      if( fabs(I) > epsabs ) {
        // printf("\nadding gaunt entry\n");
        // fflush(stdout);
 
        gaunt.push_back(I);
        gLpval.push_back(indexL(lp,mp));
      }
    }
  }

  // final gaunt index position
  gindex.push_back(gaunt.size());
  
  // free gauss-legendre table of abscissae and weights
  gsl_integration_glfixed_table_free(gltable);

}

// gsl callback function from gaunt coefficient integrand
struct gauntparam { int l1,m1,l2,m2,lp,mp; };
double gauntintegrand(double x, void *p) { 
  struct gauntparam *gp = (gauntparam*)(p);
 
  const double tpi = 2.0*M_PI;

  // generate normalized associated Legendre polynomials
  double sphL1 = gsl_sf_legendre_sphPlm(gp->l1,gp->m1,x);
  double sphL2 = gsl_sf_legendre_sphPlm(gp->l2,gp->m2,x);
  double sphLp = gsl_sf_legendre_sphPlm(gp->lp,gp->mp,x); 

  return tpi*sphL1*sphL2*sphLp;
}

// use adaptive quadrature method to perform integration
void green0_t::calc_all_gaunt_adapt() {

  // setup gsl callback function for integrand 
  gsl_function gslcb;
  gslcb.function = &gauntintegrand;

  // clear any existing gaunt table
  gaunt.clear(); gLpval.clear(); gindex.clear();

  // integration error tolerances
  const double epsabs = 1e-12, epsrel = 1.e-12;
  double I, abserr, maxerr = 0.0;

  // workspace for gsl adaptive integration
  const int wssize = 500;
  gsl_integration_workspace *wspace = 
    gsl_integration_workspace_alloc(wssize);  
  
  // go through every (L1,L2) pair
  for(int l1 = 0, L1 = 0; l1 <= maxl; l1++)
  for(int m1 = -l1; m1 <= l1; m1++, L1++)
  for(int l2 = 0, L2 = 0; l2 <= maxl; l2++)
  for(int m2 = -l2; m2 <= l2; m2++, L2++) {
  
    // set pointer into gaunt array for (L1,L2)
    gindex.push_back(gaunt.size());

    // lp must satisfy the triangle inequality
    int lp_lo = abs(l1-l2), lp_hi = l1+l2;
    for(int lp = lp_lo; lp <= lp_hi; lp++)
    for(int mp = -lp; mp <= lp; mp++) {
    
      // gaunt is zero if m values don't match
      // if( m1 + m2 != mp ) continue; // <-- correct
      if( m2 != m1 + mp ) continue; // <-- agree w/ MECCA?

      // gaunt is zero if overall parity odd 
      if( (l1+l2+lp) % 2 == 1 ) continue;

      // perform adapative integration to within 'epsabs' error
      gauntparam gp = { l1,abs(m1), l2,abs(m2), lp,abs(mp) };
      gslcb.params = &gp;
      gsl_integration_qag(&gslcb, -1.0, 1.0, epsabs, epsrel, 
        wssize, GSL_INTEG_GAUSS61, wspace, &I, &abserr);

      // minus factors since gsl harmonics only for m > 0
      double sign = 1.0;
      if( m1 >= 0 && m1%2 == 1 ) sign = -sign;
      if( m2 >= 0 && m2%2 == 1 ) sign = -sign;
      if( mp >= 0 && mp%2 == 1 ) sign = -sign;
      I = sign*I;

      // append nonzero gaunt value to list
      if( fabs(I) > epsabs ) {
        gaunt.push_back(I);
        gLpval.push_back(indexL(lp,mp));
      }

      if( abserr > maxerr ) maxerr = abserr;
    }
  }

  // final gaunt index position
  gindex.push_back(gaunt.size());
  
  // free gsl adaptive integration workspace   
  gsl_integration_workspace_free(wspace);
}




