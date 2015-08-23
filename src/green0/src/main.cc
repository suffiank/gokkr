/* written by Suffian Khan Aug 2014 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <ctime>
#include "green0.h"
#include "../../util/src/loadkvp.h"

using namespace std;

typedef complex<double> cplx;
const double tpi = 2.0*M_PI;

int main(int argc, char **argv) {

  // read in all key/value pairs from input
  map<string,string> kvp;
  bool good = safeloadkvp("input.txt",&kvp);
  if( !good ) return -1;

  try {

    // configure structure constants module 
    green0_t g0module(kvp);

    // calculate structure constants
    cplx energy( cplx(atof(kvp["energy.re"].c_str()),
      atof(kvp["energy.im"].c_str())) );
  
    vec3 kpoint; 
    kpoint.x = atof(kvp["kx"].c_str());
    kpoint.y = atof(kvp["ky"].c_str());
    kpoint.z = atof(kvp["kz"].c_str());
    vector<vec3> kpoints(1,kpoint);
    g0module.setkpoints(kpoints);
  
    int matdim = g0module.get_g0dim();
    vector<cplx> Gmat(matdim*matdim);
    clock_t start, finish;
    start = clock();
    for(int i = 0; i < 1; i++)
      g0module.calc_g0mat(Gmat.data(), energy, 0);
    finish = clock();
    double elapsed = (double(finish-start)/CLOCKS_PER_SEC);
    printf("Total G0 time: %10.6f sec\n",elapsed);

    g0module.export_gmat("gmat.out",Gmat.data());
  
  } catch(string message) {
    printf(message.c_str());
    printf("\n");
    return -1;
  }

  return 0;
}
