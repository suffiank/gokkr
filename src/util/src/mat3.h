#ifndef _MAT3_H
#define _MAT3_H

#include "vec3.h"

struct mat3 {
  vec3 row[3];
  inline vec3& operator[](int i) { return row[i]; }
  inline const vec3& operator[](int i) const { return row[i]; }
};

inline mat3 operator-(const mat3& m) {
  mat3 w; w[0] = -m[0]; w[1] = -m[1]; w[2] = -m[2];
  return w;
}

inline vec3 operator*(const mat3& m, const vec3& v) {
  return vec3( dot(m.row[0],v), dot(m.row[1],v), dot(m.row[2],v) );
}

inline mat3 operator*(const mat3& a, const mat3& b) {
  mat3 w;
  for(int i = 0; i < 3; i++)
  for(int j = 0; j < 3; j++)
    w[i][j] = a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j];
  return w;
}

inline mat3 transpose(const mat3& a) {
  mat3 w;
  for(int i = 0; i < 3; i++)
  for(int j = 0; j < 3; j++) 
    w[i][j] = a[j][i];
  return w;
}

inline bool operator==(const mat3& a, const mat3& b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

inline bool operator!=(const mat3& a, const mat3& b) {
  return !(a == b);
}

#endif

