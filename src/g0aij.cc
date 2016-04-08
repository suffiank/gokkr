/* written by Suffian Khan Aug 2014 */

#include "green0.h"
using namespace std;

void green0_t::gen_aijlist() {

  // clear current aijlist list
  aijlist.clear();

  // iterate through every (i,j) site pair
  for(int i = 0; i < natom; i++)
  for(int j = 0; j < natom; j++) {

    // fix aij = basis_i - basis_j;
    vec3 aij = basis[i]+(-1.0)*basis[j];

    // if aij already exists in aijlist, pin (i,j) to it
    bool unique = true;
    for(int k = 0; k < aijlist.size(); k++)
      if( aij == aijlist[k] ) { 
        aij_sources[k].push_back(pair<int,int>(i,j));
        unique = false; break; 
      }        
    
    // if aij does not already exist, add it to aijlist
    if(unique) {
      aijlist.push_back(aij);
      vector<pair<int,int> > t(1);
      t[0].first = i; t[0].second = j;
      aij_sources.push_back(t);
    }
  }
  export_aij();
}


