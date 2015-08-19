#include "timing.h"
#include <ctime>
#include <algorithm>
using namespace std;

const int pad = 65;

void timer_t::begin(string label) {
  if(opentimers.size() > 0 ) printf("\n");
  printf("%-*s",pad,label.c_str());
  same_line = true;
  opentimers.push( pair<string,double>(label,double(clock())) );
}

void timer_t::end() {
  clock_t finish = clock();
  pair<string,double> current = opentimers.top();
  string label = current.first;
  clock_t start = current.second;
  double elapsed = (double(finish-start)/CLOCKS_PER_SEC);
  opentimers.pop();
  if(!same_line) {
    int mpad = max(0,pad-int(label.length())-6);
    printf("End '%s'%*s",label.c_str(),mpad,"");
  }
  printf("%10.6f sec\n",elapsed);
  same_line = false;
}
