#include "crystal.h"
#include <cstdlib>
using namespace std;

static mat3 idmat;
bool symmetry_t::init_std_symm = false;

static const dble tol = 1.e-06;
// static const dble pi  = M_PI;

symmetry_t::symmetry_t() {

  // load rotations once for all objects
  // warning: this may not be thread safe
  if( !init_std_symm ) {

    // identity matrix
    idmat[0] = vec3(1.0,0.0,0.0);
    idmat[1] = vec3(0.0,1.0,0.0);
    idmat[2] = vec3(0.0,0.0,1.0);

    // standard hexagonal and cubic symmetries
    fill_cube_isometries(cube_symmetry);
    fill_hex_isometries(hex_symmetry);
    init_std_symm = true;
  }
}

void symmetry_t::configure(std::map<std::string,std::string>& kvp) {

  find();
}

void symmetry_t::find_bravais_symmetry() {

  mat3 ainv = transpose(crystal->bvec);
  basegroup = "hex";

  beginning:
  int nrot; subgroup.clear();
  if( basegroup == "cube" ) nrot = 48;
  if( basegroup == "hex"  ) nrot = 24;

  int n = 0; mat3 rot;
  for(int i = 0; i < nrot; i++) {
      if( basegroup == "cube" ) rot = cube_symmetry[i];
      if( basegroup == "hex" )  rot = hex_symmetry[i];

      // rewrite rotation in lattice coordinates 
      // i.e. R' = A^-1 R A for A = basis vectors
      rot = rot * crystal->avec;
      rot = ainv * rot;

      // check whether operation preserves Bravais lattice
      // i.e. are all entries integer
      for(int j = 0; j < 3; j++) 
      for(int k = 0; k < 3; k++) {
        dble x = fmod(rot[j][k],1.0);
        if( x < 0.0 ) x += 1.0;
        if( tol < x && x < 1.0-tol )
          goto nextsymm;
      }

      // passed check, save symmetry
      subgroup.push_back(i); n++;

      nextsymm: continue;
  }

  // determine crystal system by operations survived
  system = "";
  if( basegroup == "hex" ) {
    if( n == 24 ) system = "hexagonal";
    if( n == 12 ) system = "trigonal";
    if( system == "" ) 
      { basegroup = "cube"; goto beginning; }
  }
  else {
    if( n == 48 ) system = "cubic";
    if( n == 16 ) system = "tetragonal";
    if( n ==  8 ) system = "orthorhombic";
    if( n ==  4 ) system = "monoclinic";
    if( n ==  2 ) system = "triclinic";
    if( system == "" )
      { printf("symmetry: failed to identify crystal system\n"); exit(1); }
  }
  suborder = subgroup.size();
}

void symmetry_t::find_crystal_symmetry() {

    // determine crystal system and base group
    find_bravais_symmetry();
    printf("system = %s\n",system.c_str());

    // rewrite basis in lattice coordinates in 1st cell
    mat3 ainv = transpose(crystal->bvec);

    // for each bravais symmetry found
    nsymm = 0; mat3 rot, rotl;
    for(int i = 0; i < suborder; i++) {

      if( basegroup == "cube" )
        rot = cube_symmetry[ subgroup[i] ];
      if( basegroup == "hex" ) 
        rot = hex_symmetry[ subgroup[i] ];

      // rewrite rotation in lattice coordinates 
      // i.e. R' = A^-1 R A for A = basis vectors
      rotl = rot * crystal->avec;
      rotl = ainv * rotl;

      // for each basis site
      bool is_crystal_symm = true;
      for(int j = 0; j < crystal->basis.size(); j++) {

        // apply rotation to basis site
        vec3 vec = rotl * crystal->basisl[j];
        modulo_coord(vec);
       
        // does it coincide with some equivalent site?
        bool found_match = false;
        for(int k = 0; k < crystal->basis.size(); k++)
        if( crystal->siteid[j] == crystal->siteid[k] ) {
          vec3 dvec = vec - crystal->basisl[k];
          modulo_coord(dvec);
          if( (dvec[0] < tol || tol-1.0 < dvec[0]) &&
              (dvec[1] < tol || tol-1.0 < dvec[1]) &&
              (dvec[2] < tol || tol-1.0 < dvec[2]) ) 
          { found_match = true; break; }
        }
       
        // if not, this is not a symmetry operation 
        if( !found_match )
          { is_crystal_symm = false; break; }

      }
      if( !is_crystal_symm ) continue;

      // record valid symmetry operation
      crystal_symmetry.push_back(rot); nsymm++;
    }

    // debug check: print symmetries found 
    /* const vector<mat3>& symm = crystal_symmetry;
    printf("crystal symmetries\n");
    for(int i = 0; i < nsymm; i++) {
      printf("operation = %d\n",i+1);
      printf("%20.15f %20.15f %20.15f\n",
        symm[i][0][0],symm[i][0][1],symm[i][0][2]);
      printf("%20.15f %20.15f %20.15f\n",
        symm[i][1][0],symm[i][1][1],symm[i][1][2]);
      printf("%20.15f %20.15f %20.15f\n",
        symm[i][2][0],symm[i][2][1],symm[i][2][2]);
    } */
}

void symmetry_t::fill_hex_isometries(vector<mat3>& rotmat) {

    // zero entries be default
    rotmat.resize(24);

    // identity matrix
    int n = 0;
    rotmat[n++] = idmat;

    // rotations of 60 deg about z-axis
    vec3 axis = vec3(0.0, 0.0, 1.0);
    fill_rotation_axis(rotmat[n], pi/3.0, axis);
    fill_period(rotmat, n); n++;
 
    // rotations of 180 deg about in-plane axes
    axis = vec3(1.0, 0.0, 0.0);
    fill_rotation_axis(rotmat[n++], pi, axis);
  
    axis = vec3(cos(pi/6.0), sin(pi/6.0), 0.0);
    fill_rotation_axis(rotmat[n++], pi, axis);

    axis = vec3(cos(pi/3.0), sin(pi/3.0), 0.0);
    fill_rotation_axis(rotmat[n++], pi, axis);

    axis = vec3(0.0, 1.0, 0.0);
    fill_rotation_axis(rotmat[n++], pi, axis);

    axis = vec3(-cos(pi/3.0), sin(pi/3.0), 0.0);
    fill_rotation_axis(rotmat[n++], pi, axis);

    axis = vec3(-cos(pi/6.0), sin(pi/6.0), 0.0);
    fill_rotation_axis(rotmat[n++], pi, axis);

    // debug check: print all operations 
    /* printf("hexagonal symmetries\n");
    for(int i = 0; i < 12; i++) {
      printf("operation = %d\n",i+1);
      printf("%20.15f %20.15f %20.15f\n",
        rotmat[i][0][0],rotmat[i][0][1],rotmat[i][0][2]);
      printf("%20.15f %20.15f %20.15f\n",
        rotmat[i][1][0],rotmat[i][1][1],rotmat[i][1][2]);
      printf("%20.15f %20.15f %20.15f\n",
        rotmat[i][2][0],rotmat[i][2][1],rotmat[i][2][2]);
    } */

    // apply inversion to previous operations
    for(int i = 0; i < 12; i++)
      rotmat[n++] = -rotmat[i];    
}

void symmetry_t::fill_cube_isometries(vector<mat3>& rotmat) {

    // there are 48 cubic symmetries
    rotmat.resize(48);

    // identity matrix
    int n = 0;
    rotmat[n++] = idmat;

    // rotations about face-center to face-center
    vec3 axis = vec3( 1.0, 0.0, 0.0 );
    fill_rotation_axis(rotmat[n], pi/2.0, axis);
    fill_period(rotmat, n); n++;

    axis = vec3( 0.0, 1.0, 0.0 );
    fill_rotation_axis(rotmat[n], pi/2.0, axis);
    fill_period(rotmat, n); n++;

    axis = vec3( 0.0, 0.0, 1.0 );
    fill_rotation_axis(rotmat[n], pi/2.0, axis);
    fill_period(rotmat, n); n++;

    // rotations about edge-center to edge-center
    axis = vec3( 1.0, 1.0, 0.0 );
    fill_rotation_axis(rotmat[n++], pi, axis);
 
    axis = vec3( -1.0, 1.0, 0.0 );
    fill_rotation_axis(rotmat[n++], pi, axis);

    axis = vec3( 1.0, 0.0, 1.0 );
    fill_rotation_axis(rotmat[n++], pi, axis);

    axis = vec3( 1.0, 0.0, -1.0 );
    fill_rotation_axis(rotmat[n++], pi, axis);
    
    axis = vec3( 0.0, 1.0, 1.0 );
    fill_rotation_axis(rotmat[n++], pi, axis);

    axis = vec3( 0.0, 1.0, -1.0 );
    fill_rotation_axis(rotmat[n++], pi, axis);
 
    // rotations about body-diagonal 
    axis = vec3( 1.0, 1.0, 1.0 );
    fill_rotation_axis(rotmat[n], 2.0*pi/3.0, axis);
    fill_period(rotmat, n); n++;

    axis = vec3( -1.0, 1.0, 1.0 );
    fill_rotation_axis(rotmat[n], 2.0*pi/3.0, axis);
    fill_period(rotmat, n); n++;

    axis = vec3( 1.0, -1.0, 1.0 );
    fill_rotation_axis(rotmat[n], 2.0*pi/3.0, axis);
    fill_period(rotmat, n); n++;

    axis = vec3( 1.0, 1.0, -1.0 );
    fill_rotation_axis(rotmat[n], 2.0*pi/3.0, axis);
    fill_period(rotmat, n); n++;

    // debug check: print all operations 
    /* printf("cubic symmetries\n");
    for(int i = 0; i < 24; i++) {
      printf("operation = %d\n",i+1);
      printf("%10.5f %10.5f %10.5f\n",
        rotmat[i][0][0],rotmat[i][0][1],rotmat[i][0][2]);
      printf("%10.5f %10.5f %10.5f\n",
        rotmat[i][1][0],rotmat[i][1][1],rotmat[i][1][2]);
      printf("%10.5f %10.5f %10.5f\n",
        rotmat[i][2][0],rotmat[i][2][1],rotmat[i][2][2]);
    } */
 
    // apply inversion to previous operations
    for(int i = 0; i < 24; i++)
      rotmat[n++] = -rotmat[i];
}

void symmetry_t::fill_rotation_axis(mat3& rotmat, dble theta, const vec3& axis) {

  // see wikipedia
  dble c = cos(theta), s = sin(theta);
  dble m = sqrt(axis.x*axis.x+axis.y*axis.y+axis.z*axis.z);
  dble x = axis.x/m, y = axis.y/m, z = axis.z/m;

  rotmat[0][0] = (1.0-x*x)*c + x*x;
  rotmat[0][1] = -z*s-x*y*c + x*y;
  rotmat[0][2] = y*s - x*z*c + x*z;

  rotmat[1][0] = z*s - x*y*c + x*y;
  rotmat[1][1] = (1.0-y*y)*c + y*y;
  rotmat[1][2] = -x*s -y*z*c + y*z;

  rotmat[2][0] = -y*s - x*z*c + x*z;
  rotmat[2][1] = x*s - y*z*c + y*z;
  rotmat[2][2] = (1.0-z*z)*c + z*z;
}
  
void symmetry_t::fill_period(vector<mat3>& rotmat, int& n) {

  int n0 = n;

  mat3 mat = rotmat[n] * rotmat[n0];
  while(mat != idmat) { 
    rotmat[++n] = mat;
    mat = rotmat[n] * rotmat[n0];
  }
}


