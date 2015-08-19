/* written by Suffian Khan August 2014 */

#include "green0.h"
#include <fstream>
using namespace std;

void green0_t::export_gmat(string filename, const cplx* Gij) {

  if(log_times) timer.begin("Export G0 matrix");

  ofstream file(filename.c_str());
  char buff[256];
  int natom = basis.size();
  int stride = natom*numL;
  for(int i = 0; i < natom; i++)
  for(int j = 0; j < natom; j++) {
 
    sprintf(buff,"# i = %2i, j = %2i",i,j);
    file << buff << endl;
    sprintf(buff,"# l1,m1  l2,m2  structure constant 'Gij'"); 
    file << buff << endl;
    for(int l1 = 0; l1 <= maxl; l1++)
    for(int m1 = -l1; m1 <= l1; m1++) { 
      for(int l2 = 0; l2 <= maxl; l2++)
      for(int m2 = -l2; m2 <= l2; m2++) {

        // warning: switch to column-major ordering
        int L12 = i*numL*stride + j*numL;
        L12 += indexL(l2,m2)*stride + indexL(l1,m1);

        sprintf(buff,"  %2i %2i  %2i %2i  %22.15E %22.15E",
          l2,m2,l1,m1, Gij[L12].real(), Gij[L12].imag());
        file << buff << endl;
      }
      file << endl;
    }
    file << endl;
  }
  file.close();

  if(log_times) timer.end();
}

void green0_t::export_lattice(string filename, const vector<vec3>& latt) {
  ofstream file(filename.c_str());
  char buff[128];
  for(int i = 0; i < latt.size(); i++) {
    sprintf(buff,"%4i  %22.15E %22.15E %22.15E",
      i,latt[i].x,latt[i].y,latt[i].z);
    file << buff << endl;
  }
  file.close();
}

void green0_t::export_aij() {

  ofstream file("aij.out");
  char buff[256];
  for(int n = 0; n < aijlist.size(); n++) {
    sprintf(buff,"%4i  %22.15E %22.15E %22.15E",
      n,aijlist[n].x,aijlist[n].y,aijlist[n].z);
    file << buff << endl; 
    sprintf(buff,"---- "); file << buff;
    for(int m = 0; m < aij_sources[n].size(); m++) {
      int i = aij_sources[n][m].first;
      int j = aij_sources[n][m].second;
      sprintf(buff,"(%i,%i) ",i,j); 
      file << buff;
    }
    file << endl;
  }
  file.close();
}

void green0_t::export_rlylm() {
 
  ofstream file("rlylm.out");
  char buff[128];
  for(int i = 0; i < aijlist.size(); i++) {
    sprintf(buff,"# begin r = %22.15E %22.15E %22.15E:",
      aijlist[i].x,aijlist[i].y,aijlist[i].z);
    file << buff << endl;
    for(int j = 0; j < rlatt.size(); j++) {
      sprintf(buff,"# r = %22.15E %22.15E %22.15E",
        aijlist[i].x,aijlist[i].y,aijlist[i].z);
      file << buff << endl;
      sprintf(buff,"# R = %22.15E %22.15E %22.15E",
        rlatt[j].x,rlatt[j].y,rlatt[j].z);
      file << buff << endl;
      sprintf(buff,"# l  m                  (r+R)^l Ylm(r+R)");
      file << buff << endl;
      for(int l = 0; l <= maxlp; l++)
      for(int m = 0; m <= l; m++) {
        sprintf(buff," %2i %2i %22.15E  %22.15E",
          l,m,rlylm_latt[i][j][indexL(l,m)].real(), 
          rlylm_latt[i][j][indexL(l,m)].imag());
        file << buff << endl;
      }
      file << endl;
    }
  }
  file.close();
}

void green0_t::export_klylm() {
 
  ofstream file("klylm.out");
  char buff[128];
  for(int i = 0; i < kpoints.size(); i++) {
    sprintf(buff,"# begin k = %22.15E %22.15E %22.15E:",
      kpoints[i].x,kpoints[i].y,kpoints[i].z);
    file << buff << endl;
    for(int j = 0; j < klatt.size(); j++) {
      sprintf(buff,"# k = %22.15E %22.15E %22.15E",
        kpoints[i].x,kpoints[i].y,kpoints[i].z);
      file << buff << endl;
      sprintf(buff,"# K = %22.15E %22.15E %22.15E",
        klatt[j].x,klatt[j].y,klatt[j].z);
      file << buff << endl;
      sprintf(buff,"# l  m                  (k+K)^l Ylm(k+K)");
      file << buff << endl;
      for(int l = 0; l <= maxlp; l++)
      for(int m = 0; m <= l; m++) {
        sprintf(buff," %2i %2i %22.15E  %22.15E",
          l,m,klylm_latt[i][j][indexL(l,m)].real(), 
          klylm_latt[i][j][indexL(l,m)].imag());
        file << buff << endl;
      }
      file << endl;
    }
  }
  file.close();
}

void green0_t::export_intfac() {

  ofstream file("intfac.out");
  char buff[256]; 
  for(int i = 0; i < aijlist.size(); i++) {
    sprintf(buff,"# begin r = %22.15E %22.15E %22.15E:",
      aijlist[i].x,aijlist[i].y,aijlist[i].z); 
    file << buff << endl; 
    for(int j = 0; j < rlatt.size(); j++) {
      sprintf(buff,"# r = %22.15E %22.15E %22.15E",
        aijlist[i].x,aijlist[i].y,aijlist[i].z); 
      file << buff << endl;
      sprintf(buff,"# R = %22.15E %22.15E %22.15E",
        rlatt[j].x,rlatt[j].y,rlatt[j].z);
      file << buff << endl;
      sprintf(buff,"# l  int[ x^2l exp(-x^2 (r+R)^2 + z/4x^2)," 
        " x = sqrt(eta)/2 .. inf"); 
      file << buff << endl;
      for(int l = 0; l <= maxlp; l++) {
        sprintf(buff," %2i %22.15E %22.15E",
          l,intfac[i][j][l].real(),intfac[i][j][l].imag()); 
        file << buff << endl;
      }
      file << endl;
    }
  }
  file.close();
}

void green0_t::export_gaunt() {
 
  ofstream file("gaunt.out");
  char buff[256];
  sprintf(buff,"# l1,m1  l2,m2  lp,mp    gaunt");
  file << buff << endl;
  for(int Lj = 0; Lj < numL; Lj++)
  for(int Li = 0; Li < numL; Li++) {
    
    int L12 = numL*Li + Lj;
    int ib = gindex[L12], ie = gindex[L12+1];
    for(int i = ib; i != ie; i++) {
      sprintf(buff,"  %2i %2i  %2i %2i  %2i %2i  %22.15E",
        indexl(Li),indexm(Li),
        indexl(Lj),indexm(Lj),
        indexl(gLpval[i]),indexm(gLpval[i]),gaunt[i]);
      file << buff << endl;      
    }
  }
  file.close();
}


