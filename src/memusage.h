// written by Suffian Khan
// circa 3/23/2016

#ifndef _MEMSAMPLE_H
#define _MEMSAMPLE_H

#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <utility>
#include <cstdio>

struct memseg_t {
  const char* filename;
  int fileline;
  double timestamp;
  size_t size;
  void *ploc;
};

struct memclass_t {
  std::string label;
  std::list<memseg_t> seglist;
  std::unordered_map<void*,memseg_t> seghash;
};

struct memcnt_t {
  int cnt; long long size;
};

struct memusage_t {

  memusage_t(); 
  ~memusage_t();

  // wrappers to track malloc and free
  void* malloc(size_t size, std::string label, 
    const char* file, int line); 
  void* calloc(size_t num, size_t size, std::string label, 
    const char* file, int line);
  void* realloc(void* ploc, size_t size, std::string label, 
    const char* file, int line);
  void free(void* ploc, const char* file, int line);

  // query total VM from /proc/self/status
  int getVM();
  void log(); 
  void dump(); 

  // run as background thread
  void autolog(double interval);

private:

  void update_snapshot(void *pold, void *pnew, size_t size, 
    std::string label, const char* file, int line);
  void error(std::string msg);
  
  std::vector<memclass_t> snapshot;
  std::unordered_map<std::string,memcnt_t> counter;
  bool do_counting; int count_limit;

  std::FILE *tracef, *logf;
  double start;
  int nprinted;
  pthread_t tid;
};

extern memusage_t memusage;

#ifdef OVERRIDE_MALLOC

#define malloc(s)    memusage.malloc((s),"tracked",__FILE__,__LINE__)
#define calloc(n,s)  memusage.calloc((n),(s),"tracked",__FILE__,__LINE__)
#define realloc(p,s) memusage.realloc((p),(s),"tracked",__FILE__,__LINE__)
#define free(p)      memusage.free((p),__FILE__,__LINE__)

#endif

#endif
