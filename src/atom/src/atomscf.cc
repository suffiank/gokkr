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

  // set number of grid points
  int Z;
  if( argc == 1 )
    Z = 1;
  else
    Z = atoi(argv[1]);

  // set which atom to solve
  atom_t atom;
  atom.Z = Z;

  // fill an initial set of energy levels
  atom.set_logarithmic_grid(10000, 0.0001, 200.0);
  atom.fill_madelung_orbitals();
  atom.fill_thomas_fermi_potential();
  atom.solve_self_consistent();
}
