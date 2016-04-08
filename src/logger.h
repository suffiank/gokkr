#ifndef _LOGGER_H
#define _LOGGER_H

#include <string>
#include <algorithm>
#include <set>
#include <cstdlib>
#include <cstdarg>

struct log_t {

  void focus(const std::string& category) { 
    active = false;
    if( filter.find(category) != filter.end() )
      active = true;
  }

  void logf(const char *fmt, ...) {
    if( !active ) return;

    va_list args; va_start(args, fmt);
    std::vprintf(fmt, args);
    va_end(args);
  }
  
  void vlogf(const char *fmt, va_list args) {
    if( !active ) return;
    std::vprintf(fmt, args);
  }
  
  void enable(const std::string& category)
    { filter.insert(category); }

  void disable(const std::string& category)
    { filter.erase(category); }

  log_t() : active(false) {}

private:
  std::set< std::string > filter;
  bool active;
};

extern log_t logger;

// abbreviation to above
inline void logf(const char *fmt, ...) 
  { va_list args; va_start(args, fmt); logger.vlogf(fmt, args); va_end(args); }


#endif
