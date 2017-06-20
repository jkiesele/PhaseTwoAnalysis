//////////////////////////////////////////////////
//
//  StdArg.hpp: main's command line parser
//
//  Created: 04/29/2003, Andriy Zatserklyaniy
//  
/////////////////////////////////////////////////

#ifndef StdArg_hpp
#define StdArg_hpp

#include <map>
#include <sstream>
#include <iostream>

using std::cout;    using std::cerr;    using std::endl;
using std::string;  using std::map;

class StdArg {
//
//  Intended to parse command line arguments
//  in form of key-value pairs and flags in any order like:
//
//           flag   key-value     key-value      flag
//  ./a.out -debug -in input.dat -out output.dat -mc
//
//  See example of usage below.
//
public:
  struct BadInput  // exception handler
  {
    BadInput(const string& mess, const string& key) {
      cerr<< "*** ERROR StdArg: " << mess << key <<endl;
    }
    BadInput(const string& mess, const string& key, const string& value) {
      cerr<< "*** ERROR StdArg: " << mess << ". Key-value pair: " << key << " " << value <<endl;
    }
  };
  
  struct FlagMap: public std::map<string, int> {
    FlagMap& operator << (const string& flag) {
      insert( std::pair<string, int>(flag,0) );
      return *this;
    }
  };
  struct KeyValueMap: public std::map<string, const char*> {
    KeyValueMap& operator << (const string& key) {
      insert( std::pair<string, const char*>(key,0) );
      return *this;
    }
  };
  
  FlagMap&     flags;  // user interface for input flags
  KeyValueMap& keys;   // user interface for input keys
private:
  int _argc;
  char** _argv;
  // input, output map for flags ('1' means 'has only one field: key' )
  FlagMap     _i1, _o1;
  // input, output map for pairs key-value ('2' means 'has two fields: key and value' )
  KeyValueMap _i2, _o2;
  bool _check_minus;
public:
  StdArg(int argc, char** argv, bool check_minus=true):
    flags(_i1), keys(_i2), _argc(argc), _argv(argv), _check_minus(check_minus)
  {}
  void Process() {
    int iarg = 0;
    while (++iarg < _argc) {
      string key(_argv[iarg]);
      // first look in defined flags
      if (_i1.find(key) != _i1.end())
      { // found
        if (_o1[key] != 0) throw BadInput("defined twice flag: ", key);
        else {
          _o1[key]++;
          continue;   // next argument
        }
      }
      // not found in flags, look in defined keys for key-value
      if (_i2.find(key) != _i2.end())
      { // found
        if (_o2.find(key) != _o2.end()) throw BadInput("defined twice key: ", key);
        if (iarg+1 < _argc) _o2[key] = _argv[++iarg];  // store value and point to next key
        else throw BadInput("no value for key ", key);
      }
      else throw BadInput("no such key: ", key);
    }
  }
  void ShowFlags() const {
    // flags were read from input
    for (FlagMap::const_iterator p = _o1.begin(); p!=_o1.end(); p++) {
      cout<< p->first <<endl;
    }
  }
  void ShowKeys() const {
    // key-value pairs were read from input
    for (KeyValueMap::const_iterator p = _o2.begin(); p!=_o2.end(); p++) {
      cout<< p->first <<" "<< p->second <<endl;
    }
  }
  // flag for flags
  bool Flag(const string& flag) const { return _o1.find(flag) != _o1.end(); }
  // key for key-value pair
  bool Key(const string& key) const {
    return _o2.find(key) != _o2.end();
  }
  // doesn't throw exception
  const char* Value(const string& key) {
    if (!Key(key)) return 0;
    else           return _o2[key];
  }

  // template convertor for (almost) arbitrary type
  //
  // usage: Type n = instance.Get<Type>("-key");
  //
  template<class T> T Get(const string& key) {
    T val;
    char ch;
    if (!Key(key)) throw BadInput("no such key: ", key);
    std::istringstream buf(_o2[key]);
    if (!(buf >> val))  throw BadInput("conversation error", key, Value(key));
    else if (buf >> ch) throw BadInput("conversation error", key, Value(key));  // ok if trailing '\0'
    return val;
  }
};
//
// Specializations: should be defined out of class declaration
//
template<>  // specialization for const char*
const char* StdArg::Get<const char*>(const string& key) {
  if (!Key(key)) throw BadInput("no such key: ", key);
  if (_check_minus && *Value(key) == '-')
    throw BadInput("value should not start from '-'", key, Value(key));
  return Value(key);
}
template<>  // specialization for string: allows few words in string
string StdArg::Get<string>(const string& key) {
  if (!Key(key)) throw BadInput("no such key: ", key);
  if (_check_minus && *Value(key) == '-')
    throw BadInput("value should not start from '-'", key, Value(key));
  return string(Value(key));
}
template<>  // specialization for bool
bool StdArg::Get<bool>(const string& key) {
  bool val;
  char ch;
  if (!Key(key)) throw BadInput("no such key: ", key);
  if (strcmp(_o2[key],"true") ==0) return true;
  if (strcmp(_o2[key],"false")==0) return false;
  std::istringstream buf(_o2[key]);
  if (!(buf >> val))  throw BadInput("conversation error", key, Value(key));
  else if (buf >> ch) throw BadInput("conversation error", key, Value(key));  // ok if trailing '\0'
  return val;
}


#endif  // StdArg_hpp
