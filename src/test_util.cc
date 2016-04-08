#include <iostream>
#include <vector>
using namespace std;

#include "loadkvp.h"
#include "logger.h"
#include "mathfn.h"

int main(int argc, char **argv) {

  map<string,string> kvp;
  bool good = safeloadkvp("input.txt",&kvp);
  if( !good ) return -1;

  int n = 11;
  vector<dble> x(n), w(n);
  fill_gauss_legendre_table(n, x.data(), w.data());

  logger.enable("stdout");
  logger.focus("stdout");
  logf("%3s %20s%20s\n","i","x","w");
  for(int i = 0; i < n; i++)
    logf("%3i %20.15f%20.15f\n",i+1,x[i],w[i]);
}
