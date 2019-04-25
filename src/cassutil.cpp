/** 
 * LoxcodeR
 * Cassette utils
 */

#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <sstream>
#include "pack.h"

using namespace Rcpp;

std::vector<int> get_cass_vec_single(std::string c){
  std::stringstream s(c);
  std::vector<int> out;
  int k = 0;
  while(!s.eof()){
    s >> k;
    out.push_back(k);
  }
  return out;
}

bool is_valid_single(std::vector<int> c){
  if(c.size() != 13 &&
     c.size() != 9 && 
     c.size() != 7 &&
     c.size() != 5 &&
     c.size() != 3){
    // non-allowed size
    return false;
  }
  std::set<int> odd; std::set<int> even;
  for(int i = 0; i < c.size(); i++){
    if(i % 2 == 0){
      if(abs(c[i]) % 2 != 1) return false; // odd elt, must have odd values
      odd.insert(abs(c[i])); 
    }else{
      if(abs(c[i]) % 2 != 0) return false; // even elt, must have even values
      even.insert(abs(c[i])); 
    }
  }
  if(odd.size() != (c.size()/2 + 1) || even.size() != c.size()/2) return false;
  return true;
}

//' @export
// [[Rcpp::export]]
std::vector<std::vector<int> > get_cass_vec(std::vector<std::string> c){
  std::vector<std::vector<int> > out(c.size());
  for(int i = 0; i < c.size(); i++) out[i] = get_cass_vec_single(c[i]);
  return out;
}

template<typename T>
std::vector<bool> is_valid_impl(T c){
  return std::vector<bool>(); // DO NOT USE
}

template<>
std::vector<bool> is_valid_impl(std::vector<std::vector<int> > c){
  std::vector<bool> out(c.size());
  for(int i = 0; i < c.size(); i++) out[i] = is_valid_single(c[i]);
  return out;
}

template<>
std::vector<bool> is_valid_impl(std::vector<std::string> c){
  std::vector<std::vector<int> > c_int = get_cass_vec(c);
  return is_valid_impl(c_int);
}


//' @export
// [[Rcpp::export]]
std::vector<bool> is_valid(SEXP c){
  switch(TYPEOF(c)){
    case INTSXP:
    case REALSXP:
      return is_valid_impl<std::vector<std::vector<int> > >(as<std::vector<std::vector<int> > >(c));
      break;
    case STRSXP:
      return is_valid_impl<std::vector<std::string> >(as<std::vector<std::string> >(c));
      break;
  }
  return std::vector<bool>();
}
