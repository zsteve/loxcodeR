/**
 * LoxcodeR
 * Cassette distances
 */

#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <sstream>
#include <tuple>
#include "pack.h"
#include "cassutil.h"

using namespace Rcpp;

int size_to_size_idx[] = {-1, -1, -1, 4, -1, 3, -1, 2, -1, 1, -1, -1, -1, 0};

// [[Rcpp::export]]
void wrapper_fill_tables(){
  // wrapper
  fill_tables();
  cerr << __FUNCTION__ << ": tables filled" << endl;
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<int> > get_oes(std::vector<int>& k){
  std::vector<int> o(k.size()/2 + 1);
  std::vector<int> e(k.size()/2);
  std::vector<int> s(k.size());
  for(int i = 0; i < k.size(); i++){
     if(i % 2 == 0) o[i/2] = (k[i]-1)/2; // odd elt
     else e[i/2] = (k[i]/2)-1;

     s[i] = k[i] < 0 ? 1 : 0;
  }
  return std::make_tuple(o, e, s);
}

template<typename T>
std::vector<long long> pack_impl(T c){ return std::vector<long long>(); /* DO NOT USE */ }

template<>
std::vector<long long> pack_impl(std::vector<std::vector<int> > c){
  std::vector<long long> out(c.size());
  for(int i = 0; i < c.size(); ++i){
    int size_idx = -1;
    if(c[i].size() > 13) size_idx = -1;
    else size_idx = size_to_size_idx[c[i].size()];
    if(size_idx != -1){    
      std::tuple<std::vector<int>, std::vector<int>, std::vector<int> > p = get_oes(c[i]);
      long long packed = pack(pack_odd(&(std::get<0>(p))[0], size_idx),
                          pack_even(&(std::get<1>(p))[0], size_idx),
                          pack_sign(&(std::get<2>(p))[0], size_idx), size_idx);
      out[i] = packed;
    }else{
      out[i] = -1;
    }

  }
  return out;
}

template<>
std::vector<long long> pack_impl(std::vector<std::string> c){
  std::vector<std::vector<int> > c_int = get_cass_vec(c);
  return pack_impl(c_int);
}

//' @export
// [[Rcpp::export]]
std::vector<long long> pack(SEXP c){
  switch(TYPEOF(c)){
    case VECSXP:
      return pack_impl<std::vector<std::vector<int> > >(as<std::vector<std::vector<int> > >(c));
      break;
    case STRSXP:
      return pack_impl<std::vector<std::string> >(as<std::vector<std::string> >(c));
      break;
  }
  return std::vector<long long>();
}
