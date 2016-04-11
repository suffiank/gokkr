// written by Suffian Khan
// circa 3/23/2016

#ifndef _TIMING_H
#define _TIMING_H

#include <ctime>
#include <string>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <pthread.h>

// timed segment representation
struct timeseg_t { 
 
  std::string label;
  int count; bool open;
  double accum, clockedat;
};

struct mytimer_t {

  mytimer_t() : tid( pthread_t(-1) ) {} 
  ~mytimer_t() { if( tid != pthread_t(-1) ) pthread_cancel(tid); }

  // dummy functions
  void begin(std::string label) {}
  void end(std::string label) {}
  void end() {}

  void dump();
  void autolog(double interval);

  // hash table of timers
  static const int hashsize = 104729;
  timeseg_t seghash[hashsize];
  int numtimers;
  pthread_t tid;
};

extern mytimer_t timer;

#ifdef PROFILE

// hash algorithm from LoL engine
// see http://lolengine.net/blog/2011/12/20/cpp-constant-string-hash
#define H1(s,i,x)   (x*65599u+(uint8_t)s[(i)<std::strlen(s)?std::strlen(s)-1-(i):std::strlen(s)])
#define H4(s,i,x)   H1(s,i,H1(s,i+1,H1(s,i+2,H1(s,i+3,x))))
#define H16(s,i,x)  H4(s,i,H4(s,i+4,H4(s,i+8,H4(s,i+12,x))))
#define H64(s,i,x)  H16(s,i,H16(s,i+16,H16(s,i+32,H16(s,i+48,x))))
#define H256(s,i,x) H64(s,i,H64(s,i+64,H64(s,i+128,H64(s,i+192,x))))

#define HASH(s)    ((uint32_t)(H256(s,0,0)^(H256(s,0,0)>>16)))

#define begin_timer(_label) { \
  const int id = HASH(_label) % timer.hashsize; \
  if( timer.seghash[id].count == 0 ) \
    timer.seghash[id].label = _label; \
  timer.seghash[id].count++; \
  timer.seghash[id].open = true; \
  timer.seghash[id].clockedat = double(clock()); \
}

#define end_timer(_label) { \
  const int id = HASH(_label) % timer.hashsize; \
  const double elapse = \
    (double(clock())-timer.seghash[id].clockedat)/CLOCKS_PER_SEC; \
  timer.seghash[id].open = false; \
  timer.seghash[id].accum += elapse; \
  timer.seghash[id].clockedat = elapse; \
} 

#else

#define begin_timer(label) {}
#define end_timer(label) {}

#endif


#endif
