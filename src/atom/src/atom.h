#ifndef _ATOM_H_
#define _ATOM_H_

#include <vector>
#include <map>
#include <string>
#include "../../util/src/mathfn.h"

struct orbital_t {
  int n;   // principle quantum number
  int l;   // angular momentum 
  int occ; // occupancy (i.e. filling, not degeneracy)
  dble E;  // energy

  xfunc wavef; // (radial wave function)/radius
  xfunc sigma; // (radial charge density)
};

struct atom_t {
  int Z; // atomic number
  int N; // number of grid positions

  dble mu, xi, xf, dx; // logarthmic grid parameters
  xfunc r;     // radial grid positions
  xfunc rV;    // radius * (electron potential energy)
  xfunc sigma; // radial charge density (i.e. 4 pi r^2 rho)

  int nlvl;
  std::vector<orbital_t> orbital;

  bool verbose;

  // common interface
  atom_t() : verbose(false) {}
  atom_t(std::map<std::string,std::string>& kvp) { configure(kvp); }
  void configure(std::map<std::string,std::string>& kvp);

  // methods
  void set_logarithmic_grid(int N, dble rmin, dble rmax);

  void fill_madelung_orbitals();
  void fill_thomas_fermi_potential();
  void fill_coloumb_potential();
  void fill_square_well(dble U0, dble R0);
  void fill_harmonic_oscillator(dble omega);

  int  recompute_orbital(int lvl);
  void solve_orbitals();
  void sort_orbitals();

  void recompute_charge();
  void recompute_potential();
  int  solve_self_consistent();

  // regular & irregular solution according to Zeller def.
  void solve_reg_fn(int l, dble rad, cplx E, zfunc *Rl, cplx *tl, dble *dl);
  void solve_irr_fn(int l, dble rad, cplx E, zfunc *Sl);
};

#endif
