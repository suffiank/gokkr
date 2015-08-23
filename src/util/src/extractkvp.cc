#include "extractkvp.h"
#include <cstdlib>
#include <cstdio>
#include <vector>

#define USE_MATH_DEFINES
#include <cmath>

using namespace std;

const double tpi = 2.0*M_PI;

void extractkvp(map<string,string>& kvp, const string& key, 
  double &value, double defaultval) {

  // safe extraction of real value from key/value list

  char *e;
  if( kvp.find(key) != kvp.end() ) {
    value = strtod(kvp[key].c_str(), &e);
    if(*e != '\0')
      throw string("error: unable to interpret '"+key+
        "'; '"+kvp[key]+"' invalid float");
  }
  else if( defaultval != require_key_dbl )
    value = defaultval;
  else
    throw string("error: unable to find required key '"+key+"'");
}

void extractkvp(map<string,string>& kvp, const string& key, 
  int &value, int defaultval) {

  // safe extraction of integer value from key/value list

  char *e;
  if( kvp.find(key) != kvp.end() ) {
    value = strtol(kvp[key].c_str(),&e,0);
    if(*e != '\0')
      throw string("error: unable to interpret key '"+key+
        "'; '"+kvp[key]+"' invalid integer");
  }
  else if( defaultval != require_key_int )
    value = defaultval;
  else
    throw string("error: unable to find required key '"+key+"'");
}

void extract_lattvec(map<string,string>& kvp, mat3& avec) {

  /* the Bravais vectors should be specified in the input as:
   
     crystal.lattice:
       x=A1X y=A1Y z=A1Z
       x=A2X y=A2Y z=A2Z
       x=A3X y=A3Y z=A3Z

     this results in key/value pairs:

       crystal.lattice.x=A1X
       crystal.lattice.x2=A2X
       crystal.lattice.x3=A3X
       etc.
  */


  try {
    extractkvp(kvp, "crystal.lattice.x", avec[0].x, require_key_dbl);
    extractkvp(kvp, "crystal.lattice.y", avec[0].y, require_key_dbl);
    extractkvp(kvp, "crystal.lattice.z", avec[0].z, require_key_dbl);
    // avec[0] = tpi*avec[0];
  
    extractkvp(kvp, "crystal.lattice.x2", avec[1].x, require_key_dbl);
    extractkvp(kvp, "crystal.lattice.y2", avec[1].y, require_key_dbl);
    extractkvp(kvp, "crystal.lattice.z2", avec[1].z, require_key_dbl);
    // avec[1] = tpi*avec[1];
  
    extractkvp(kvp, "crystal.lattice.x3", avec[2].x, require_key_dbl);
    extractkvp(kvp, "crystal.lattice.y3", avec[2].y, require_key_dbl);
    extractkvp(kvp, "crystal.lattice.z3", avec[2].z, require_key_dbl);
    // avec[2] = tpi*avec[2];
  }
  catch(string msgin) {
    throw string("error: unable to read Bravais vectors\n"+msgin);
  }
}

void extract_basis(map<string,string>& kvp, vector<vec3>& basis) {

  /* the basis vectors should be specified in input file as:
   
     crystal.basis:
       x=X1 y=Y1 z=Z1
       x=X2 y=Y2 z=Z2
       x=X3 y=Y3 z=Z3
       ...
       x=XN y=YN z=ZN

     this results in key/value pairs:

       crystal.basis.x=X1
       crystal.basis.x2=X2
       crystal.basis.x3=X3
       ...
       crystal.basis.xN=XN.

     and similarly for y and z
  */ 

  // clear previous basis
  basis.clear();

  // find the number of atoms by searching for largest
  // key of the form 'basis.x'+N
  // note: the below search process is slow, O(N log N)

  int n = 1; char buff[32]; 
  string header = "crystal.basis.", cnt;
  do{
    n++; sprintf(buff,"x%i",n); cnt = buff;
  } while( kvp.find(header+cnt) != kvp.end() );

  if(n-1 == 1 && kvp.find(header+"x") == kvp.end() )
   throw string("error: no basis sites defined");

  int natom = n-1;

  // confirm 'basis.y'+N and 'basis.z'+N match

  n = 1;
  do{
    n++; sprintf(buff,"y%i",n); cnt = buff;
  } while( kvp.find(header+cnt) != kvp.end() );

  if(n-1 != natom)
    throw string("error: number of basis sites is ill-defined");

  n = 1;
  do{
    n++; sprintf(buff,"z%i",n); cnt = buff;
  } while( kvp.find(header+cnt) != kvp.end() );

  if(n-1 != natom)
    throw string("error: number of basis sites is ill-defined");

  // assign basis vectors

  basis.resize(natom);
  string ax = "x", ay = "y", az = "z";
  for(n = 0; n < natom; n++) {
    sprintf(buff,"");
    if(n > 0) sprintf(buff,"%i",n+1);
    cnt = buff;

    try {
      extractkvp(kvp, header+ax+cnt, basis[n].x, require_key_dbl);
      extractkvp(kvp, header+ay+cnt, basis[n].y, require_key_dbl);
      extractkvp(kvp, header+az+cnt, basis[n].z, require_key_dbl);
      // basis[n] = tpi*basis[n];
    }
    catch(string msgin) {
      throw string("error: unable to read basis sites\n"+msgin);
    }
  }
}


