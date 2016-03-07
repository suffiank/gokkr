#ifndef _YLM_H
#define _YLM_H

#include "vec3.h"
#include "mathfn.h"

void calc_ylm_norm(const int maxl, double *Alm);
void calc_plm(const int maxl, double *Plm, const double& cth, const double& sth);
void calc_vlylm(const int maxl, const double *Alm, cplx *vlylm, const vec3& v);

#endif


