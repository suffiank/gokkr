#include "timer.h"
#include <algorithm>
#include <ctime>
using namespace std;

timer_t timer;
typedef map< string, pair<double,double> >::iterator iter_t;

void timer_t::begin(string label) {

  iter_t T = times.find(label);
  if( T != times.end() )
    T->second.second = double(clock());
  else
    // warning: this line is broken! leads to seg fault
    T->second = pair<double,double>(0.0, double(clock()));
  lastlabel = label;
}

void timer_t::end(string label) {

  clock_t finish = clock();
  iter_t T = times.find(label);
  if( T != times.end() ) {
    clock_t start = T->second.second;
    double elapsed = (double(finish-start)/CLOCKS_PER_SEC);
    T->second.first += elapsed; T->second.second = elapsed;
  }
  else
    throw string("timer_t:: attempt to end() non-existant clock");
}

void timer_t::clear(string label) {

  iter_t T = times.find(label);
  if( T != times.end() ) times.erase(T);
}

double timer_t::last(string label) {
  
  iter_t T = times.find(label);
  if( T != times.end() )
    return T->second.second;
  else
    throw string("timer_t:: attempt to last() non-existant clock");
}

double timer_t::total(string label) {

  iter_t T = times.find(label);
  if( T != times.end() )
    return T->second.first;
  else
    throw string("timer_t:: attempt to total() non-existant clock");
}
