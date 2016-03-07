/* written by Suffian Khan Aug 2014 */

#ifndef _GREEN0_HH
#define _GREEN0_HH

#include <string>
#include <vector>
#include <map>
#include <complex>
#include <ctime>
#include "../../util/src/mat3.h"
#include "indexlm.h"
#include "../../util/src/timer.h"

typedef std::complex<double> cplx;
class crystal_t;

struct green0_t {

  green0_t() : crystal(NULL) {}
  green0_t(std::map<std::string,std::string>& kvp);
  void configure(std::map<std::string,std::string>& kvp);
  void bind(const crystal_t* cryst);

  // cache options, note: enable multiple modes with bit-wise OR
  typedef int cache_t;
  static const cache_t CACHE_NONE=0x00, CACHE_RSPACE=0x01, CACHE_KSPACE=0x02;

  // calculation of k-space G0: Fourier transform or Ewald method
  typedef int calc_t;
  static const calc_t USE_FOURIER=1, USE_EWALD=2;

  void setmaxl(const int maxl);
  void setewaldparam(const double eta, const double maxR, const double maxK);
  void setcrystal(const mat3& avec, const std::vector<vec3>& basis);
  void setcachemode(const cache_t& cachemode);
  void setcalcmode(const calc_t& calcmode);
  void setkpoints(const std::vector<vec3>& kpoints);

  int get_g0dim();
  void calc_g0mat(cplx *g0mat, const cplx& z, const vec3& k);
  void calc_g0mat(cplx *g0mat, const cplx& z, const int kindex);
  void export_gmat(std::string filename, const cplx *g0mat);

private:

  // for timing calculations
  bool log_times;

  // crystal description
  // todo: replace w/ crystal*
  mat3 avec, bvec;
  std::vector<vec3> basis;
  int natom; double vol;
  const crystal_t* crystal;

  // size of basis expansion
  int maxl, numL;

  // ewald parameters
  int maxlp, numLp; // Lp means L' from Ewald sum
  double eta, maxR, maxK;

  // real and reciprocal lattice vectors
  std::vector<vec3> rlatt, klatt;
  void genlattice(std::vector<vec3> *latt, const mat3& avec, 
    const double& rad);
  void export_lattice(std::string filename, const std::vector<vec3>& latt);

  // calculate (i,j) block of green0 in k-space
  void calc_g0ij_kspace(cplx *g0ij, const int stride, const cplx& z, 
    const int aij, const int kindex);
  calc_t calcmode;

  // remember last energy requested to avoid some recomputation
  cplx lastz;

  // Dlm r- and k-space lattice sums
  // note: different versions exist based on the caching scheme

  void calc_dlm_fourier(cplx *Dlm, const cplx& z, const cplx& p,
    const vec3& r, const vec3& k); 
  void calc_dlm_fourier(cplx *Dlm, const cplx& z, const cplx& p,
    const int aij, const vec3& k); 

  void add_dlm_ewald_rsum(cplx *Dlm, const cplx& z, const cplx& p,
    const vec3& r, const vec3& k);
  void add_dlm_ewald_rsum(cplx *Dlm, const cplx& z, const cplx& p,
    const int aij, const vec3& k);
  void add_dlm_ewald_ksum(cplx *Dlm, const cplx& z, const cplx& p,
    const vec3& r, const vec3& k);
  void add_dlm_ewald_ksum(cplx *Dlm, const cplx& z, const cplx& p,
    const vec3& r, int kindex);

  cplx calc_d00(const cplx& z);
  cplx d00;

  // data structures for related to 
  // harmonic polynomials (R+r)^l Y_lm(R+r) [or (K+k)^l Y_lm(K+k)]:

  // represents 'r' points for which (R+r)^l Y_lm(R+r) are saved
  // it is always the case r = aij = basis_i - basis_j
  std::vector<vec3> aijlist;
  void gen_aijlist();
  void export_aij();;
  cache_t cachemode;

  // aij_sources[n][*] gives (i,j) pairs with aij = aijlist[n]
  std::vector< std::vector< std::pair<int,int> > > aij_sources;

  // represents (R+r)^l Y_lm(R+r)
  // first index is r, second index R, third index L = lm
  // indices correspond to lists in 'aijlist' and 'rlatt'
  std::vector< std::vector< std::vector<cplx> > > rlylm_latt;

  // represents 'k' points for which (K+k)^l Y_lm(K+k) are saved
  std::vector<vec3> kpoints;

  // represents (K+k)^l Y_lm(K+k)
  // first index is k, second index K, third index L = lm
  // indices correspond to lists in 'kpoints' and 'klatt'
  std::vector< std::vector< std::vector<cplx> > > klylm_latt;

  // represents hankel function h+_l[p(r+R)] for Fourier transform
  // first index is r=aij, second index R, third index l
  std::vector< std::vector< std::vector<cplx> > > hankel;
  void calc_hankel(cplx *hl, const cplx& x);
  void calc_all_hankel(const cplx& z);
  void resize_all_hankel();

  // calculation of spherical harmonics:
  // Alm = sqrt( (2l+1)/4pi (l-m)!/(l+m)! ) 
  // note: plm storage index increments l before m
  std::vector<double> Alm;

  void calc_all_rlylm();
  void export_rlylm();

  void calc_all_klylm();
  void export_klylm();

  // integration factor int[x^2l exp(-x^2 (R+r)^2 + b^2/4z^2) dx]
  // first index r, second index R, third index l
  std::vector< std::vector< std::vector<cplx> > > intfac;
  void calc_intfac(cplx *intfac, const cplx& z, const double& r);
  void resize_all_intfac();
  void calc_all_intfac(const cplx& z);
  void export_intfac();

  // represents gaunt coefficients int[ Y_L1(t) Y_L2(t) Y_Lp(t)* sin(t) dt]
  // note: linear storage scheme for efficiency
  std::vector<double> gaunt; // gaunt values in order (L1,L2,L')
  std::vector<int> gLpval;   // the L' corresponding to gaunt[i]
  std::vector<int> gindex;   // start index of (L1,L2) in gaunt array
  void calc_all_gaunt();
  void calc_all_gaunt_fixed();
  void calc_all_gaunt_adapt();
  void export_gaunt();
};

#endif
