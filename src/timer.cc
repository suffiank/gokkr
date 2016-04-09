// written by Suffian Khan
// circa 3/23/2016

#include "timer.h"
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <pthread.h>
#include <unistd.h>
using namespace std;

mytimer_t timer;

// strip profiling in production code
#ifndef PROFILE

void mytimer_t::dump() {}

#else

void mytimer_t::dump() {

  FILE* file = fopen("timing.out","w");
    
  int maxlabel = 10, maxtime = 12;
  int maxcount = 8, maxid = 8;
  for(int i = 0; i < hashsize; i++)
    if( seghash[i].label.length() != 0 ) {

    char buff[48]; int len;

    if( seghash[i].label.length() > maxlabel ) 
      maxlabel = seghash[i].label.length();

    snprintf(buff, 48, "%i", i); 
    len = strlen(buff);
    if( len > maxid ) maxid = len;
 
    snprintf(buff, 48, "%i", seghash[i].count); 
    len = strlen(buff);
    if( len > maxcount ) maxcount = len;
    
    snprintf(buff, 48, "%.5f", seghash[i].accum); 
    len = strlen(buff);
    if( len > maxtime ) maxtime = len;
  }

  fprintf(file,"#%-*s%*s%*s%*s%*s\n", 
    maxlabel, " label", ++maxid, "hash", ++maxcount, "count", 
    ++maxtime, "total(s)", maxtime, "average(s)");

  double curr = double(clock());
  for(int i = 0; i < hashsize; i++)
    if( seghash[i].label.length() != 0 ) {
      timeseg_t& ts = seghash[i];
      
      fprintf(file,"%-*s%*d%*d%*.5f%*.5f\n", 
        maxlabel+1, ts.label.c_str(), maxid, i,
        maxcount, ts.count, maxtime, 
        ts.accum + (ts.open?(curr-ts.clockedat)/CLOCKS_PER_SEC:0), 
        maxtime, ts.accum/ts.count);
    }
 
  fclose(file);
}

void* autolog_timer_main(void* arg) {

  double interval = *( (double*)arg );
  int microsec = int(1000000.*interval);
  while(true) 
   { timer.dump(); usleep(microsec); }
  return NULL;
} 

void mytimer_t::autolog(double interval) {

  static double arg = interval;
  if( tid != pthread_t(-1) ) pthread_cancel(tid);
  pthread_create(&tid, NULL, autolog_timer_main, (void*)(&arg));
} 

#endif
