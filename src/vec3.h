#ifndef _VEC3_H
#define _VEC3_H

#include <algorithm>
#include <cmath>
#include "mathfn.h"

#define _vec3tol_ 1e-13;

struct vec3 { 
  dble x,y,z; 

  vec3() : x(0.0), y(0.0), z(0.0) {}
  vec3(dble _x, dble _y, dble _z) : x(_x), y(_y), z(_z) {}
  inline dble& operator[](int i) { 
    switch(i) { 
      case 0: return x;
      case 1: return y;
      case 2: return z;
    }; throw;
  }
  inline const dble& operator[](int i) const { 
    switch(i) { 
      case 0: return x;
      case 1: return y;
      case 2: return z;
    }; throw;
  } 
  static const dble tol;
};


inline dble dot(const vec3& u, const vec3& v) {
  return u.x*v.x + u.y*v.y + u.z*v.z;
}
inline vec3 cross(const vec3& u, const vec3& v) {
  return vec3( 
    u.y*v.z - u.z*v.y,
    u.z*v.x - u.x*v.z,
    u.x*v.y - u.y*v.x );
} 
inline dble mag(const vec3 v) {
  return sqrt(dot(v,v));
}
inline vec3 operator-(const vec3& v) {
  return vec3(-v.x,-v.y,-v.z);
}
inline vec3 operator+(const vec3& u, const vec3& v) {
  return vec3(u.x+v.x, u.y+v.y, u.z+v.z);
}
inline vec3 operator-(const vec3& u, const vec3& v) {
  return vec3(u.x-v.x, u.y-v.y, u.z-v.z);
}
inline vec3 operator*(const dble a, const vec3& v) {
  return vec3(a*v.x, a*v.y, a*v.z);
}
inline vec3 operator/(const vec3& v, const dble a) {
  return vec3(v.x/a, v.y/a, v.z/a);
}
// note: detect equality based on fractional difference
inline bool operator==(const vec3& u, const vec3& v) {
  dble mag = 0.0;
  mag = std::max(mag, std::fabs(u.x));
  mag = std::max(mag, std::fabs(u.y));
  mag = std::max(mag, std::fabs(u.z));
  mag = std::max(mag, std::fabs(v.x));
  mag = std::max(mag, std::fabs(v.y));
  mag = std::max(mag, std::fabs(v.z));
  const dble tol = _vec3tol_;
  if( mag == 0.0 ) return true;
  if( std::fabs(u.x-v.x)/mag < tol && 
      std::fabs(u.y-v.y)/mag < tol && 
      std::fabs(u.z-v.z)/mag < tol )
    return true;
  return false;
}
struct sortvec3_mag_less {
  bool operator()(const vec3& u, const vec3& v) {
    return mag(u) < mag(v);
  }
};
struct sortvec3_xyz_less {
  bool operator()(const vec3& u, const vec3& v) {
    const dble tol = _vec3tol_;
    if(u.x < v.x - tol) return true;
    if(u.x > v.x + tol) return false;
    if(u.y < v.y - tol) return true;
    if(u.y > v.y + tol) return false;
    if(u.z < v.z - tol) return true;
    if(u.z > v.z + tol) return false;
    return false;
  }
};

inline void modulo_coord(vec3& point) {

  for(int i = 0; i < 3; i++) {
    point[i] = fmod(point[i],1.0);
    if( point[i] < 0.0 ) point[i] += 1.0;
  }
}

#undef _vec3tol_
#endif
