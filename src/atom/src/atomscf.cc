/* 
  A version of the Herman & Skillman program for 
  non-relativistic Hartree-Fock-Slater solutions of a single atom 
    by Suffian Khan, started: 04/24/2011, working: 05/10/2015

  exclusively using atomic Rydberg units
  (i.e. hbar = 2m = e^2/2 = a0 = 1, c = 2/alpha)
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
using namespace std;

#include "../../util/src/mathfn.h"
#include "../../util/src/loadkvp.h"
#include "atom.h"

char speclett(int l) {
  switch(l) {
    case 0: return 's';
    case 1: return 'p';
    case 2: return 'd';
    case 3: return 'f';
    default: 
      return char(int('g')+l-4);
  };
}

int main(int argc, char *argv[]) {

  // read in all key/value pairs from input
  map<string,string> kvp;
  bool good = safeloadkvp("input.txt",&kvp);
  if( !good ) return -1;

  // set which atom to solve
  atom_t atom(kvp);

  printf("  solving Z = %d\n",atom.Z);
  printf("grid points = %d\n",atom.N);
  printf(" radial min = %-20.15f\n",atom.r[0]);
  printf(" radial max = %-20.15f\n",atom.r[atom.N-1]);

  // fill an initial set of energy levels
  atom.fill_madelung_orbitals();
  atom.fill_thomas_fermi_potential();
  atom.solve_self_consistent();
}
