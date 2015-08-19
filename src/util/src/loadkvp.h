/* written by Suffian Khan Aug 2014 */

#ifndef _LOADKVP_H
#define _LOADKVP_H

#include <fstream>
#include <string>
#include <cstring>
#include <cstdio>
#include <map>
#include <vector>

struct loadkvp {

  struct readerror {
    std::string message;
    std::string badline;
    int badpos;
  };
  
  loadkvp(std::string filename, std::map<std::string,std::string>* kvp);

private:

  struct readstate {
    std::string file;
    std::string line;
    int linec, pos;
  } rs;
  
  /* error handling */
  int msglen; char* message;

  /* do not allow instantiation of object */
  loadkvp() {}
  loadkvp(const loadkvp& lkvp) {}

  void splitkeyvalue(std::string word, std::string *key, std::string *value);
  std::string grabword(std::string line, int *pos, int *nws = 0);
  void throwerror(int pos = -1);

};

#endif
