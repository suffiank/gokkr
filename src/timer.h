#ifndef _TIMING_H
#define _TIMING_H

#include <ctime>
#include <cstdio>
#include <string>
#include <map>
#include <utility>

struct timer_t {
  void begin(std::string label);
  void end(std::string label);
  void end() { end(lastlabel); }
  void clear(std::string label);
  double last(std::string label);
  double total(std::string label);

private:
  // first is cumulative time
  // second is last timed segment (if closed)
  //   or begin clock (if open)
  std::map< std::string, std::pair<double,double> > times;
  std::string lastlabel;
};

extern timer_t timer;

#endif
