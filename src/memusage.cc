// written by Suffian Khan
// circa 3/23/2016

#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <pthread.h>
#include <unordered_map>
#include <mutex>
#include <time.h>
#include <dlfcn.h>
#include "memusage.h"
using namespace std;

static const bool use_hash = true;
memusage_t memusage;

// strip profiling in production code
#ifndef PROFILE

memusage_t::memusage_t() {}
memusage_t::~memusage_t() {}

int memusage_t::getVM() { return 0; }
void memusage_t::log() {} 
void memusage_t::autolog(double interval) {} 
void memusage_t::error(string msg) {}

void* memusage_t::malloc(size_t size, string label, 
  const char* file, int line) { return NULL; } 
void* memusage_t::calloc(size_t num, size_t size, string label, 
  const char* file, int line) { return NULL; }
void* memusage_t::realloc(void *ploc, size_t size, string label, 
  const char* file, int line) { return NULL; }
void memusage_t::free(void *ploc, const char* file, int line) {}

#else

#ifdef OVERRIDE_MALLOC
#undef malloc
#undef calloc
#undef realloc
#undef free
#endif

void memusage_t::update_snapshot(void *pold, void *pnew, size_t size, 
  string label, const char* file, int line) {

  // grab some statistics
  double current= double(clock());
  double timestamp = (current-start)/CLOCKS_PER_SEC;
  memseg_t block = {file, line, timestamp, size, pnew};
  string oldkey; size_t oldsize = -1; 
  
  // hash on file/line combo to check counters
  memcnt_t* p_cnt; int cnt = 0; 
  if( do_counting ) {

    string key = string(file) + ":" + to_string(line);
    p_cnt = &counter[key]; cnt = ++(p_cnt->cnt);
    if( cnt > count_limit ) {
      if( cnt == count_limit+1 )
        fprintf(logf,"%14s<< warning: suppressing %21s >>\n",
          "",key.c_str());
    }
  }
  
  // search and remove existing block
  int seek = 0;
  vector<memclass_t>::iterator mc;
  if( pold != NULL ) {

    // iterate over user defined memory class
    for(mc = snapshot.begin(); mc != snapshot.end(); mc++) {
 
      if( use_hash ) {
        // --hash based search   
        unordered_map<void*,memseg_t>& table = mc->seghash;
        unordered_map<void*,memseg_t>::iterator ms;
        ms = table.find(pold);
 
        // erase memory block if found
        if( ms != table.end() ) { 
          oldkey = string(ms->second.filename)+":"+
            to_string(ms->second.fileline);
          oldsize = ms->second.size; 
          table.erase(ms); break; 
        }
      }
      else {
        // --list based search   
        list<memseg_t>::iterator ms, H, T;
        list<memseg_t>& blocks = mc->seglist;
  
        // search alternating from front and back 
        ms = H = blocks.begin(); T = blocks.end();
        bool onhead = true;
        while( H != T && ms->ploc != pold ) 
         { ms = (onhead? --T:++H); onhead = !onhead; seek++; }
   
        // erase memory block if found
        if( ms != blocks.end() ) {
          oldkey = string(ms->filename)+":"+
            to_string(ms->fileline);
          oldsize = ms->size; 
          blocks.erase(ms); break; 
        }
      }
    }

    // report possible leak if block unfound
    if( mc == snapshot.end() )
      fprintf(logf,"%11.1f%15s%6i%10s%10d%15p\n", 
        timestamp, file, line, "lost?", 0, pold); 
  }
  
  // check to add block to a user defined class
  if( pnew != NULL ) {

    // search for user defined memory class
    for(mc = snapshot.begin(); mc != snapshot.end(); mc++)
      if( mc->label == label ) break;
  
    // make new class if it doesn't exist
    if( mc == snapshot.end() ) {
      
      memclass_t newclass;
      newclass.label = label;
      if( use_hash ) 
        newclass.seghash.insert(pair<void*,memseg_t>(block.ploc,block));
      else
        newclass.seglist.push_front(block);

      snapshot.push_back(newclass);
    }

    // otherwise just insert this block
    else 
      if( use_hash ) 
        mc->seghash.insert( pair<void*,memseg_t>(block.ploc,block) );
      else
        mc->seglist.push_front(block);
  }

  // update aggregate size and check suppression
  if( do_counting ) {
    if( pold != NULL ) counter[oldkey].size -= (long long)oldsize;
    if( pnew != NULL ) p_cnt->size += (long long)size;
    if( cnt > count_limit ) return;
  }

  // log memory events
  if( pold != NULL )
    fprintf(logf,"%11.1f%15s%6i%10s%10d%15p %-8d\n", 
      timestamp, file, line, "destroy", (int)oldsize, pold, seek); 

  if( pnew != NULL )
    fprintf(logf,"%11.1f%15s%6i%10s%10d%15p %-8d\n", 
      timestamp, file, line, "create", (int)size, pnew, cnt);

  fflush(logf);
}


void* memusage_t::malloc(size_t size, string label, 
  const char* file, int line) {
  
  void* ploc = ::malloc(size);
  update_snapshot(NULL, ploc, size, label, file, line);
  return ploc;
}

void* memusage_t::calloc(size_t num, size_t size, string label, 
  const char* file, int line) {

  void* ploc = ::calloc(num, size);
  update_snapshot(NULL, ploc, num*size, label, file, line);
  return ploc;
}

void* memusage_t::realloc(void *ploc, size_t size, string label, 
  const char* file, int line) { 

  void *pnew = ::realloc(ploc, size);
  update_snapshot(ploc, pnew, size, label, file, line);
  return pnew;
}

void memusage_t::free(void *ploc, const char* file, int line) {

  ::free(ploc);
  update_snapshot(ploc, NULL, 0, "", file, line); 
}

memusage_t::memusage_t() {

  // open files and print headers
  tracef = fopen("memtrace.out","w");
  fprintf(tracef, "#%17s%18s\n","cpu-time(s)","virtual-mem(kb)");

  logf = fopen("memevent.out","w");
  fprintf(logf,"#%10s%15s%6s%10s%10s%15s %-8s\n",
    "cpu-time","file","line","event","size","address","etc");

  // initialization
  start = double(clock());
  nprinted = 0; tid = -1;
  count_limit = 100;
  do_counting = true;
}

memusage_t::~memusage_t() {
  if( tid != -1 ) pthread_cancel(tid);
  fclose(tracef); fclose(logf);
}

void memusage_t::error(string msg) {
  printf("%s\n",msg.c_str()); exit(-1);
}

int memusage_t::getVM() { 

  FILE *stat = fopen("/proc/self/status","r");
  char line[256];
  int vm = 0;
  
  bool found = false;
  while( fgets(line,256,stat) ) 
    if( strncmp(line,"VmSize:",7) == 0) {

      found = true;
      int l = 0, r = strlen(line)-3;
      while( !isdigit(line[l]) ) l++;
      line[r] = '\0';  
      vm = atoi( &line[l] );
    }

  if(!found) error("memusage:: failed to sample /proc/self/status");
 
  fclose(stat);
  return vm;
}

void memusage_t::dump() {

  FILE* file = fopen("memsnap.out","w");
  fprintf(file, "#%11s%10s%21s\n","size","invk","file_line"); 

  long long total = 0;
  unordered_map<string,memcnt_t>::iterator it;
  for(it = counter.begin(); it != counter.end(); it++) {
    if( it->second.size > 0 ) 
      fprintf(file,"%12lld%10d%21s\n",
        it->second.size, it->second.cnt, it->first.c_str());
    total += it->second.size;
  }
  fprintf(file,"\n# total = %lld bytes\n",total);

  fclose(file);
}

void memusage_t::log() {

  static mutex critsec;  
  critsec.lock();

  // output any new class labels
  vector<memclass_t>::iterator mc;
  if( nprinted < snapshot.size() ) {
    fprintf(tracef, "#%17s%18s","cpu-time(s)","virtual-mem(by)");
    for(mc = snapshot.begin(); mc != snapshot.end(); mc++) 
      fprintf(tracef, "%18s", mc->label.c_str()); 
    fprintf(tracef,"\n");
    nprinted = snapshot.size();
  }
 
  // append current consumption
  double current = double(clock());
  double elapsed = (current-start)/CLOCKS_PER_SEC;
  fprintf(tracef,"%18.5f%18lld",elapsed,1024*((long long)getVM()));

  for(mc = snapshot.begin(); mc != snapshot.end(); mc++) {
    
    long long total = 0;
    if( use_hash ) { 
      unordered_map<void*,memseg_t>& blocks = mc->seghash;
      unordered_map<void*,memseg_t>::iterator ms;
      for(ms = blocks.begin(); ms != blocks.end(); ms++)
        total += (long long)ms->second.size;
    }
    else {
      list<memseg_t>& blocks = mc->seglist;
      list<memseg_t>::iterator ms;
      for(ms = blocks.begin(); ms != blocks.end(); ms++)
        total += (long long)ms->size;
    }
    fprintf(tracef,"%18lld",total); 
  }
  fprintf(tracef,"\n");

  fflush(tracef);
  critsec.unlock();
}

void* autolog_main(void* arg) {

  double interval = *( (double*)arg );
  int microsec = int(1000000.*interval);
  while(true) 
   { memusage.log(); memusage.dump(); usleep(microsec); }
  return NULL;
} 

void memusage_t::autolog(double interval) {

  static double arg = interval;
  if( tid != -1 ) pthread_cancel(tid);
  pthread_create(&tid, NULL, autolog_main, (void*)(&arg));
} 

#endif
