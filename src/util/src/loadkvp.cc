/* written by Suffian Khan Aug 2014 */

#include "loadkvp.h"
#include <fstream>
#include <iostream>
using namespace std;

loadkvp::loadkvp(string filename, map<string,string>* kvp) {

  // maintain list of sections 
  std::vector<std::string> sections;
  bool haskey;
 
  // make room for potential error message 
  msglen = 256;
  message = new char[msglen];

  // reset list of sections 
  sections.clear(); haskey = false;
   
  // for each line from input 
  rs.file = filename;
  ifstream file(filename.c_str()); 
  string line; rs.linec = 0;
  while( getline(file, line) ) {
    rs.line = line; rs.linec++; 

    // get the first word on line 
    int pos = 0, nws, indent = 0;
    string word = grabword(line,&pos,&indent);
    rs.pos = indent;
    
    // ignore comment lines 
    if( word.length() == 0 || word[0] == '#' )
      continue;

    // check indentation is correct
    if( indent%2 != 0 || indent/2 > sections.size() ) {
      sprintf(message, "indentation is incorrect");
      throwerror();
    }

    // remove any sections based on indent 
    if( indent/2 < sections.size() ) {

      // make sure we assigned some key to this section 
      if(!haskey) {
        sprintf(message,"section '%s' has no keys",
          sections[sections.size()-1].c_str());
        throwerror(); 
      }

      // remove sections 
      int nremove = sections.size() - indent/2;
      for(int j = 0; j < nremove; j++)
        sections.pop_back();
    }

    // if first word ends with a ':'
    if( word[word.length()-1] == ':' ) {

      // strip ':' from word 
      word = word.substr(0,word.length()-1);

      // check that the section name exists 
      if(word == "") { 
        sprintf(message,"failed to read section name");
        throwerror(); 
      }

      // add name to section list 
      sections.push_back(word);
      haskey = false;

      // grab next word 
      word = grabword(line,&pos);
      rs.pos = pos-word.length();
    }

    // otherwise for each space delimited word on line 
    while(word != "") {

      // skip remaining comment 
      if(word[0] == '#') break;

      // identify key/value pair/
      string key, value;
      splitkeyvalue(word, &key, &value);

      // prepend section headers 
      string keyf = "";
      for(int j = 0; j < sections.size(); j++)
        keyf += sections[j] + ".";
      keyf += key;

      // do not add keys that already exist 
      /* if( kvp->find(keyf) != kvp->end() ) {
        sprintf(message,"key '%s' already exists",keyf.c_str());
        throwerror();
      } */
      // if the key already exists
      if( kvp->find(keyf) != kvp->end() ) {
        
        // find the last increment
        int j = 1; char buff[32]; string cnt;
        do{
          j++; sprintf(buff,"%i",j); cnt = buff;
        } while( kvp->find(keyf+cnt) != kvp->end() );
        
        // add increment to key string
        keyf = keyf + cnt;         
      }

      // add to dictionary 
      (*kvp)[keyf] = value;
      haskey = true;
     
      // move to next word
      word = grabword(line, &pos);
      rs.pos = pos-word.length(); 
    }
  }

  delete [] message;
}


// 'word' should be of the form 'key=value' 
void loadkvp::splitkeyvalue(string word, string *key, string *value) {

  int wordl = word.length();

  // find location of '=' to split key/value pair
  int i = 0;
  while(i < wordl && word[i] != '=') i++;
  if(i == wordl) {
    sprintf(message,"failed to interpret '%s' as key/value pair",
      word.c_str());
    throwerror();
  }
 
  // check key begins with alphabetic character
  if( wordl > 0 && !isalpha(word[0]) ) {
    sprintf(message,"key '%s' does not begin with an alphabetic "
      "char", word.substr(0,i).c_str());
    throwerror();
  }

  // split key and value strings
  *key   = word.substr(0,i);
  *value = word.substr(i+1,wordl-i-1);
  
  // check value is not empty 
  if( value->length() == 0 ) {
    sprintf(message,"key '%s' has no value definition",
      word.substr(0,i).c_str());
    throwerror(rs.pos+i);
  }
}

// grab next word from 'line' as delimited by spaces
// single quotes may be used to retain spaces 
// on entry, 'pos' indicates where to begin search
// on exit,  'pos' indicates the character after word
// 'nws' indicates the number of whitespace skipped 
string loadkvp::grabword(string line, int *pos, int *nws) {

  int p; string word = "";

  // find first not whitespace 
  for(p = *pos; p < line.length() && line[p] == ' '; p++);
  if(nws) *nws = p-*pos;

  // find next whitespace 
  bool quote = false, ws = false;
  while( (!ws || quote) && p < line.length() ) {
    switch( line[p] ) {
      case ' ': 
        if(quote) word.append(1,line[p]);
        ws = true; break;
      case '\'': 
        quote = !quote; 
        ws = false; break;
      default:  
        word.append(1,line[p]);
        ws = false; break;
    }
    p++;
  }

  // set position and return word 
  *pos = (p < line.length()? p-1:p);
  return word;
}

void loadkvp::throwerror(int pos) {

  if( pos < 0 ) 
    pos = rs.pos;

  char buffer[256];
  sprintf(buffer,"%s:%i:%i: error: ", rs.file.c_str(), rs.linec, pos+1);
  strcat(buffer,message);

  readerror re = {buffer, rs.line, pos};
  throw re;
}

bool safeloadkvp(string filename, map<string,string>* kvp_in) {

  map<string,string>& kvp = *kvp_in;

  // read in all key/value pairs from input
  try {
    loadkvp(filename,&kvp);

    // if requested to export key/value pairs
    if( kvp.find("export.kvp") != kvp.end() &&
        kvp["export.kvp"] == "true" ) {
       
      // print all key/value pairs to 'kvp.out'
      ofstream file("kvp.out"); char buff[256];
      sprintf(buff,"%30s %30s\n", "key", "value"); file << buff;
      sprintf(buff,"%30s %30s\n", "---", "-----"); file << buff;
      map<string,string>::iterator it = kvp.begin();
      for(; it != kvp.end(); it++) {
        sprintf(buff,"%30s %30s\n", it->first.c_str(), it->second.c_str());
        file << buff;
      } 
      file.close();
    }
  }  
  catch(loadkvp::readerror re) {
    cout << re.message << endl; 
    if(re.badline != "") cout << re.badline << endl;
    if(re.badpos >= 0) {
      for(int i = 0; i < re.badpos; i++)
        cout << ' ';
      cout << "^\n";
    }
    return false; 
  }

  return true;
}


