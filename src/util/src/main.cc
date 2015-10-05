#include <iostream>
using namespace std;

#include "loadkvp.h"
#include "logger.h"

int main(int argc, char **argv) {

  map<string,string> kvp;
  bool good = safeloadkvp("input.txt",&kvp);
  if( !good ) return -1;

  int i = 100;
  logger.enable("stdout");
  logger.focus("stdout");
  logf("My int i=%i\n",i);
}
