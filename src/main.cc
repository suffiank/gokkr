#include <iostream>
using namespace std;

#include "util/src/loadkvp.h"
#include "util/src/mat3.h"
#include "crystal/src/crystal.h"

int main(int argc, char **argv) {

  // read in all key/value pairs from input
  map<string,string> kvp;
  try {
    loadkvp("input.txt",&kvp);

    // if requested to export key/value pairs
    if( kvp.find("export.kvp") != kvp.end() &&
        kvp["export.kvp"] == "true" ) {
       
      // print all key/value pairs to 'kvp.out'
      ofstream file("kvp.out"); char buff[256];
      sprintf(buff,"%30s %30s\n", "key", "value"); file << buff;
      sprintf(buff,"%30s %30s\n", "---", "-----"); file << buff;
      map<string,string>::iterator it = kvp.begin();
      for(; it != kvp.end(); it++) {
        sprintf(buff,"%30s %30s\n", it->first.c_str(), it->second.c_str());
        file << buff;
      } 
      file.close();
    }
  }  
  catch(loadkvp::readerror re) {
    cout << re.message << endl; 
    if(re.badline != "") cout << re.badline << endl;
    if(re.badpos >= 0) {
      for(int i = 0; i < re.badpos; i++)
        cout << ' ';
      cout << "^\n";
    }
    return -1; 
  }

  crystal_t crystal(kvp);

  /* mat3 avec;
  avec[0] = vec3(0.5,0.5,0.0); 
  avec[1] = vec3(0.0,0.5,0.5); 
  avec[2] = vec3(0.5,0.0,0.5); 
  crystal.set_bravais(avec);

  vector<vec3> basis(1, vec3(0.0,0.0,0.0) );
  vector<int> siteid(1, 0);
  crystal.set_basis(basis, siteid); */

  crystal.find_symmetry();
  crystal.find_specialk(10,10,10);

  // crystal.greenfn(E, 0.01, 0.22);
  // crystal.tau00(E);
  
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
