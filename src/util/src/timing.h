#ifndef _TIMING_H
#define _TIMING_H

#include <ctime>
#include <cstdio>
#include <string>
#include <stack>
#include <utility>

struct timer_t {
  void begin(std::string label);
  void end();
private:
  std::stack< std::pair<std::string,double> > opentimers;
  bool same_line;
};

#endif
