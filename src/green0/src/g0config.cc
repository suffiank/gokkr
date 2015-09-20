/* written by Suffian Khan Aug 2014 */

#include "../../crystal/src/crystal.h"
#include "../../util/src/extractkvp.h"
#include "../../util/src/logger.h"
#include "../../util/src/ylm.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <ctime>
#include <map>

using namespace std;

green0_t::green0_t(map<string,string>& kvp) : crystal(NULL) {
  configure(kvp);
}

void green0_t::bind(const crystal_t* cryst) {
  crystal = cryst;
  avec = crystal->avec;
  basis = crystal->basis;
  natom = crystal->basis.size();
  setcrystal(this->avec,this->basis);
}

void green0_t::configure(map<string,string>& kvp) {

  // todo: nominally only green0 parameters
  //  should be read here

  // should we log times?
  log_times = false;
  string key = "green0.timing";
  if( kvp.find(key) != kvp.end() )
   if( kvp[key] == "true" ) 
     log_times = true;
   else if( kvp[key] == "false" ) 
     log_times = false;

  // configure ewald parameters
  extractkvp(kvp,"green0.ewald.eta",eta,0.338925550376917);
  extractkvp(kvp,"green0.ewald.maxR",maxR,3.0);
  extractkvp(kvp,"green0.ewald.maxK",maxK,3.36154726279432);

  // configure cache settings
  cachemode = CACHE_RSPACE;
  key = "green0.cache.rspace";
  if( kvp.find(key) != kvp.end() )
    if( kvp[key] == "true" ) 
      cachemode = cachemode|CACHE_RSPACE;
    if( kvp[key] == "false" )
      cachemode = cachemode&(~CACHE_RSPACE);

  key = "green0.cache.kspace";
  if( kvp.find(key) != kvp.end() )
    if( kvp[key] == "true" ) 
      cachemode = cachemode|CACHE_KSPACE;
    if( kvp[key] == "false" )
      cachemode = cachemode&(~CACHE_KSPACE); 

  // configure Ewald method or Fourier transform 
  calcmode = USE_EWALD;
  key = "green0.method";
  if( kvp.find(key) != kvp.end() )
   if( kvp[key] == "ewald" ) 
     calcmode = USE_EWALD;
   else if( kvp[key] == "fourier" ) 
     calcmode = USE_FOURIER;

  // configure maximum l cutoff 
  extractkvp(kvp,"approximation.maxl",maxl,3);
  setmaxl(maxl);

  // configure crystal
  if( crystal == NULL ) {
    extract_lattvec(kvp, avec);
    extract_basis(kvp, basis);
    natom = basis.size();
  }
  else {
    avec  = crystal->avec;
    basis = crystal->basis;
    natom = crystal->basis.size();
  }
  setcrystal(this->avec,this->basis);

  // print out settings
  logger.focus("green0.settings");
  logf("Green0 R-space caching: ");
  if(cachemode & CACHE_RSPACE) 
    logf("enabled\n");
  else
    logf("disabled\n");

  logf("Green0 K-space caching: ");
  if(cachemode & CACHE_KSPACE) 
    logf("enabled\n");
  else
    logf("disabled\n");

  logf("Green0 calculation scheme: ");
  if(calcmode == USE_EWALD)
    logf("Ewald method\n");
  else
    logf("Fourier transform\n");
}

void green0_t::setmaxl(const int maxl) {

  this->maxl = maxl;
  numL = (maxl+1)*(maxl+1);

  maxlp = 2*maxl;
  numLp = (maxlp+1)*(maxlp+1);

  // calculate spherical harmonic normalization factors
  Alm.resize(numLp);
  calc_ylm_norm(maxlp, Alm.data());

  // Generate gaunt coefficients 
  calc_all_gaunt();
  // export_gaunt();

  // if enabled, cache real-space values
  if(cachemode & CACHE_RSPACE) 
    if(calcmode == USE_EWALD) { 
      resize_all_intfac();
      calc_all_rlylm();
    }
    else {
      resize_all_hankel();
      calc_all_rlylm();
    }

  // if enabled, Generate k-space harmonics
  if(cachemode & CACHE_KSPACE) 
    calc_all_klylm();
}

void green0_t::setewaldparam(const double eta, const double maxR, 
  const double maxK) {
  this->eta = eta;
  this->maxR = maxR;
  this->maxK = maxK;
 
  /* reGenerate lattice, etc. */
  setcrystal(this->avec, this->basis);
}

void green0_t::setcachemode(const cache_t& cachemode) {

  // if turning on (r+R)^l Ylm(r+R) caching
  // then reGenerate real-space harmonics
  if( ~(this->cachemode&CACHE_RSPACE) & (cachemode&CACHE_RSPACE) ) {
    if(calcmode == USE_EWALD) {
      resize_all_intfac();
      calc_all_rlylm();
    }
    else {
      resize_all_hankel();
      calc_all_rlylm();
    }
  }

  // if turning on (k+K)^l Ylm(k+K) caching
  // then reGenerate k-space harmonics
  if( ~(this->cachemode&CACHE_RSPACE) & (cachemode&CACHE_RSPACE) )
    calc_all_klylm();

  // change internal cache setting
  this->cachemode = cachemode;
}

void green0_t::setcalcmode(const calc_t& calcmode) {

  this->calcmode = calcmode;
}

void green0_t::setcrystal(const mat3& avec, const vector<vec3>& basis) {

  // Generate direct lattice
  this->avec[0] = avec[0];
  this->avec[1] = avec[1];
  this->avec[2] = avec[2];
  vol = fabs(dot(avec[0],cross(avec[1],avec[2])));
  if(log_times) timer.begin("Generate R-space lattice");
  genlattice(&rlatt,avec,maxR);
  if(log_times) timer.end();

  this->basis = basis;
  natom = basis.size();

  // debugging printout
  logger.focus("green0.settings");
  logf("R-space lattice within R = %8.5E contains %i vectors.\n",
    maxR,(int)rlatt.size());

  // Generate reciprocal lattice
  bvec[0] = tpi/vol*cross(avec[1],avec[2]);
  bvec[1] = tpi/vol*cross(avec[2],avec[0]);
  bvec[2] = tpi/vol*cross(avec[0],avec[1]);
  if(log_times) timer.begin("Generate K-space lattice");
  genlattice(&klatt,bvec,maxK);
  if(log_times) timer.end();

  // debugging printout
  logf("K-space lattice within K = %8.5E contains %i vectors.\n",
    maxK,(int)klatt.size());

  // construct list of unique aij
  if(log_times) timer.begin("Generate aij = bi-bj lattice");
  gen_aijlist();
  if(log_times) timer.end();

  export_lattice("rlatt.out",rlatt);
  export_lattice("klatt.out",klatt);
  export_lattice("aij.out",aijlist);

  // if enabled, cache real-space values
  if(cachemode & CACHE_RSPACE) 
    if(calcmode == USE_EWALD) { 
      resize_all_intfac();
      calc_all_rlylm();
    }
    else {
      resize_all_hankel();
      calc_all_rlylm();
    }

  // if enabled, Generate k-space harmonics
  if(cachemode & CACHE_KSPACE) 
    calc_all_klylm();
}

int green0_t::get_g0dim() {
  return natom*numL;
}

void green0_t::calc_all_rlylm() {

  // allocated space for (r+R)^l Y_lm(r+R)
  rlylm_latt.resize(aijlist.size());
  vector<cplx> vlylm(numLp);
  for(int i = 0; i < aijlist.size(); i++) {
    rlylm_latt[i].resize(rlatt.size());
    fill(rlylm_latt[i].begin(),rlylm_latt[i].end(),vlylm);
  } 

  if(log_times) timer.begin("Generate R-space harmonics");

  // perform (r+R)^l Y_lm(r+R) calculation
  for(int i = 0; i < aijlist.size(); i++)
  for(int j = 0; j < rlatt.size(); j++) {
    if( mag(aijlist[i]+rlatt[j]) < 1.e-13 ) continue;
    vec3 rvec = aijlist[i] + rlatt[j];
    calc_vlylm(maxlp, Alm.data(), rlylm_latt[i][j].data(), rvec);
  }
  if(log_times) timer.end();

  // export_rlylm();
}

void green0_t::calc_all_klylm() {

  // Generate (k+K)^l Ylm(k+K) for k-space
  klylm_latt.resize(kpoints.size());
  vector<cplx> vlylm(numLp);
  for(int i = 0; i < kpoints.size(); i++) {
    klylm_latt[i].resize(klatt.size());
    fill(klylm_latt[i].begin(),klylm_latt[i].end(),vlylm);
  } 

  if(log_times) timer.begin("Generate K-space harmonics");

  for(int i = 0; i < kpoints.size(); i++)
  for(int j = 0; j < klatt.size(); j++) {
    if( mag(kpoints[i]+klatt[j]) < 1.e-13 ) continue;
    vec3 kvec = kpoints[i] + klatt[j];
    calc_vlylm(maxlp, Alm.data(), klylm_latt[i][j].data(), kvec);
  }
  if(log_times) timer.end();

  // export_klylm(); 
}

void green0_t::setkpoints(const vector<vec3>& kpoints) {

  this->kpoints = kpoints;
  if(cachemode & CACHE_KSPACE) 
    if(calcmode == USE_EWALD) 
      calc_all_klylm();
}

void green0_t::resize_all_hankel() {

  // resize integration factors for (r+R) and lm
  hankel.resize(aijlist.size());
  vector<cplx> ei(maxlp+1);
  for(int i = 0; i < aijlist.size(); i++) {
    hankel[i].resize(rlatt.size());
    fill(hankel[i].begin(),hankel[i].end(),ei);
  }
  lastz = nan("");
}

void green0_t::calc_all_hankel(const cplx& z) {
 
  if(log_times) timer.begin("Generate Hankel functions");
  
  cplx p = sqrt(z);
  if(p.imag() < 0.0) p = -p;

  // warning: no bounds check on 'hankel'!
  for(int i = 0; i < aijlist.size(); i++)
  for(int j = 0; j < rlatt.size(); j++) {
    vec3 rv = aijlist[i]+rlatt[j];
    double mrv = mag(rv);
    if(mrv < 1.e-10) continue;
    calc_hankel(hankel[i][j].data(),p*mrv);
  }

  if(log_times) timer.end();
}

void green0_t::resize_all_intfac() {

  // resize integration factors for (r+R) and lm
  intfac.resize(aijlist.size());
  vector<cplx> ei(maxlp+1);
  for(int i = 0; i < aijlist.size(); i++) {
    intfac[i].resize(rlatt.size());
    fill(intfac[i].begin(),intfac[i].end(),ei);
  }
  lastz = nan("");
}

void green0_t::calc_all_intfac(const cplx& z) {

  if(log_times) timer.begin("Generate Ewald integration factor");

  // warning: no bounds check on 'intfac'!
  for(int i = 0; i < aijlist.size(); i++)
  for(int j = 0; j < rlatt.size(); j++) {
    vec3 rv = aijlist[i]+rlatt[j];
    double mrv = mag(rv);
    if(mrv < 1.e-13) continue;
    calc_intfac(intfac[i][j].data(),z,mrv);
  }

  if(log_times) timer.end();

}


