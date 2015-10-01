#include <iostream>
using namespace std;

#include "../../util/src/mat3.h"
#include "../../util/src/loadkvp.h"
#include "crystal.h"

int main(int argc, char **argv) {


  /* mat3 avec;
  avec[0] = vec3(0.5,0.5,0.0); 
  avec[1] = vec3(0.0,0.5,0.5); 
  avec[2] = vec3(0.5,0.0,0.5); 
  crystal.set_bravais(avec);

  vector<vec3> basis(1, vec3(0.0,0.0,0.0) );
  vector<int> siteid(1, 0);
  crystal.set_basis(basis, siteid); */

  map<string,string> kvp;
  bool good = safeloadkvp("input.txt",&kvp);
  if( !good ) return -1;
  
  crystal_t crystal(kvp);
  // crystal.find_symmetry();
  // crystal.find_specialk(10,10,10); 
  
  /* for(int i = 0; i < crystal.symmetry.size(); i++) {
    printf("operation = %d\n",i+1);
    printf("%20.15f %20.15f %20.15f\n",
        crystal.symmetry[i][0][0],
        crystal.symmetry[i][0][1],
        crystal.symmetry[i][0][2]);
    printf("%20.15f %20.15f %20.15f\n",
        crystal.symmetry[i][1][0],
        crystal.symmetry[i][1][1],
        crystal.symmetry[i][1][2]);
    printf("%20.15f %20.15f %20.15f\n\n",
        crystal.symmetry[i][2][0],
        crystal.symmetry[i][2][1],
        crystal.symmetry[i][2][2]);
  } */
}
