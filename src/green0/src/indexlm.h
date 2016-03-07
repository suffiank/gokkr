#ifndef _INDEXLM_H
#define _INDEXLM_H

/* conversion between L and (l,m) indexing */
inline int indexL(int l, int m) { return l*(l+1) + m; }
inline int indexl(int L) {
  int l = 0;
  while( l*l <= L ) l++;
  return l-1;
}
inline int indexm(int L) { 
  int l = indexl(L);
  return L - l*(l+1);
}

#endif
