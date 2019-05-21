/**
 * LoxcodeR
 * Cassette distances
 */

#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <sstream>
#include <tuple>
#include <fstream>
#include <cassert>
#include "pack.h"
#include "cassutil.h"

using namespace Rcpp;

// oooo global variables :(((
int size_to_size_idx[] = {-1, -1, -1, 4, -1, 3, -1, 2, -1, 1, -1, -1, -1, 0};

class distmaps{
  public:
    static ifstream* origin_files[5]; // files from the origin (organised by size_idx)
    static ifstream* pair_files[2]; // only supports 13 (size_idx = 0) and 9 (size_idx = 1) cassettes
    static bool initialised;
    static bool initialised_pair;

    static void load_origin_files(std::vector<std::string> paths){
      if(initialised){
        cerr << __FUNCTION__ << ": already initialised!" << endl;
        Rcpp::stop("Already initialised");
      }
      if(paths.size() != 5){
         cerr << __FUNCTION__ << ": paths.size() was not 5! missing some distance maps..." << endl;
         Rcpp::stop("Missing files");
      }
      for(int size_idx = 0; size_idx < paths.size(); size_idx++){
        distmaps::origin_files[size_idx] = new ifstream(paths[size_idx].c_str());
        if(!distmaps::origin_files[size_idx]->is_open()){
          cerr << __FUNCTION__ << ": WARNING: failed to open file for size_idx = " << size_idx << ", path = " << paths[size_idx] << endl;
        }
      }
      initialised = true;
    }

    static void load_pair_files(std::vector<std::string> paths){
      if(initialised_pair){
        cerr << __FUNCTION__ << ": already initialised pair!" << endl;
        Rcpp::stop("Already initialised");
      }
      if(paths.size() != 2){
        cerr << __FUNCTION__ << ": paths.size() was not 2! missing some distance maps..." << endl;
        Rcpp::stop("Missing files");
      }
      for(int size_idx = 0; size_idx <= 1; size_idx++){
        distmaps::pair_files[size_idx] = new ifstream(paths[size_idx].c_str());
        if(!distmaps::pair_files[size_idx]->is_open()){
          cerr << __FUNCTION__ << ": WARNING: failed to open file for size_idx = " << size_idx << ", path = " << paths[size_idx] << endl;
        }
      }
      initialised_pair = true;
    }

    static unsigned char read_origin(long long offset, int size_idx){
      assert(size_idx < 5 && size_idx >= 0);
      origin_files[size_idx]->seekg(offset);
      return origin_files[size_idx]->get();
    }

    static unsigned char read_pair(long long offset, int size_idx){
      assert(size_idx == 0 || size_idx == 1);
      pair_files[size_idx]->seekg(offset);
      return pair_files[size_idx]->get();
    }

    ~distmaps(){
      if(initialised){
        for(int i = 0; i < 5; i++) delete origin_files[i];
        cerr << __FUNCTION__ << ": closed distance table files" << endl;
        initialised = false;
      }
    }

  private:
    distmaps() {}
};

ifstream* distmaps::origin_files[5] = {NULL};
ifstream* distmaps::pair_files[2] = {NULL};
bool distmaps::initialised = false;
bool distmaps::initialised_pair = false;
// [[Rcpp::export]]
void load_origin_files_wrapper(std::vector<std::string> paths){ distmaps::load_origin_files(paths); }

// [[Rcpp::export]]
void load_pair_files_wrapper(std::vector<std::string> paths){ distmaps::load_pair_files(paths); }

// [[Rcpp::export]]
void wrapper_fill_tables(){
  // wrapper
  fill_tables();
  cerr << __FUNCTION__ << ": tables filled" << endl;
}

int get_size_idx(int size){
  return (size > 13 || size < 0) ? -1 : size_to_size_idx[size];
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<int> > get_oes(std::vector<int>& k){
  std::vector<int> o(k.size()/2 + 1);
  std::vector<int> e(k.size()/2);
  std::vector<int> s(k.size());
  for(int i = 0; i < k.size(); i++){
     if(i % 2 == 0) o[i/2] = (abs(k[i])-1)/2; // odd elt
     else e[i/2] = (abs(k[i])/2)-1;

     s[i] = k[i] < 0 ? 1 : 0;
  }
  return std::make_tuple(o, e, s);
}

long long pack_single(std::vector<int> c, int size_idx){
  std::tuple<std::vector<int>, std::vector<int>, std::vector<int> > p = get_oes(c);
  long long packed = pack(pack_odd(&(std::get<0>(p))[0], size_idx),
                          pack_even(&(std::get<1>(p))[0], size_idx),
                          pack_sign(&(std::get<2>(p))[0], size_idx), size_idx);
  return packed;
}

template<typename T>
std::vector<long long> pack_impl(T c, std::vector<bool> v){ return std::vector<long long>(); /* DO NOT USE */ }

template<>
std::vector<long long> pack_impl(std::vector<std::vector<int> > c, std::vector<bool> v){
  std::vector<long long> out(c.size());
  for(int i = 0; i < c.size(); ++i){
    if(v[i] == false){ out[i] = NA_INTEGER;  continue;  }
    int size_idx = get_size_idx(c[i].size());
    if(size_idx != -1){
      out[i] = pack_single(c[i], size_idx);
    }else{
      out[i] = NA_INTEGER;
    }
  }
  return out;
}

template<>
std::vector<long long> pack_impl(std::vector<std::string> c, std::vector<bool> v){
  std::vector<std::vector<int> > c_int = get_cass_vec(c);
  return pack_impl(c_int, v);
}

//' Pack cassettes into cassette ID
//'
//' @param c a list of numeric vectors, or a character vector of decoded loxcodes
//' @param v vector of bool, output from validate()
//' @export
// [[Rcpp::export]]
std::vector<long long> pack(SEXP c, std::vector<bool> v){
  switch(TYPEOF(c)){
    case VECSXP:
      return pack_impl<std::vector<std::vector<int> > >(as<std::vector<std::vector<int> > >(c), v);
      break;
    case STRSXP:
      return pack_impl<std::vector<std::string> >(as<std::vector<std::string> >(c), v);
      break;
  }
  return std::vector<long long>();
}

//' @export
// [[Rcpp::export]]
std::vector<int> retrieve_dist_origin(std::vector<long long> c, std::vector<int> sizes){
  std::vector<int> out(c.size());
  if(!distmaps::initialised){
     Rcpp::stop("Not initialised");
  }
  for(int i = 0; i < c.size(); i++){
    if(c[i] != NA_INTEGER){
      // -1 indicates invalid cassette value. If not -1 then we assume that sizes is sensible...
      // must first call is_valid and pack beforehand!
      out[i] = (int)distmaps::read_origin(c[i], get_size_idx(sizes[i]))-1; // don't forget to subtract 1
    }else{
      out[i] = NA_INTEGER; // -1 for missing values
    }
  }
  return out;
}

int sign(int x){ return x > 0 ? 1 : -1; }

std::vector<int> transform_13(std::pair<std::vector<int>, std::vector<int> > c){
  std::map<int, int> elt_map;
  for(int i = 0; i < 13; i++){
    elt_map[abs(c.first[i])] = sign(c.first[i])*(i+1);
  }
  std::vector<int> out(13);
  for(int i = 0; i < 13; i++){
    out[i] = sign(c.second[i])*elt_map[abs(c.second[i])];
  }
  return out;
}

std::vector<int> transform_9(std::pair<std::vector<int>, std::vector<int> > c){
  // transform 9-element cassettes
  // transforms first cassette to 1-9, and second cassette relative to 1-9.
  std::map<int, int> elt_map;
  std::set<int> unseen_elts_e = {2, 4, 6, 8, 10, 12};
  std::set<int> unseen_elts_o = {1, 3, 5, 7, 9, 11, 13};
  for(int i = 0; i < 9; i++){
    elt_map[i+1] = c.first[i];
    if((i+1) % 2 == 1){
      // odd element
      unseen_elts_o.erase(abs(c.first[i]));
    }else{
      // even element
      unseen_elts_e.erase(abs(c.first[i]));
    }
  }
  // by this point, unseen_elts contains the 4 elements not seen in c[i]
  for(int i = 9; i < 13; i++){
    if((i+1) % 2 == 1){
      // odd element
      elt_map[i+1] = *unseen_elts_o.begin();
      unseen_elts_o.erase(unseen_elts_o.begin());
    }else{
      // even element
      elt_map[i+1] = *unseen_elts_e.begin();
      unseen_elts_e.erase(unseen_elts_e.begin());
    }
  }
  // now flip map
  std::map<int, int> rev_map;
  for(int i = 0; i < 13; i++){
    rev_map[abs(elt_map[i+1])] = sign(elt_map[i+1])*(i+1);
  }
  vector<int> out(9);
  for(int i = 0; i < 9; i++){
    out[i] = rev_map[abs(c.second[i])]*sign(c.second[i]);
  }
  return out;
}

/*
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix retrieve_dist_pair(std::vector<std::vector<int> > c, int size){
  if(size != 13 && size != 9){
    // only support distance pairs for size 9 and 13
    cerr << __FUNCTION__ << ": pairwise distance supported only for cassettes of size 13, 9" << endl;
    Rcpp::stop("Not supported");
  }
  int size_idx = get_size_idx(size);
  Rcpp::NumericMatrix out(c.size(), c.size());
  for(int i = 0; i < c.size(); i++){
    for(int j = i; j < c.size(); j++){
      if(size == 13){
        std::vector<int> x = transform_13(std::make_pair(c[i], c[j]));
        long long x_packed = pack_single(x, 0);
        int dist = distmaps::read_pair(x_packed, 0)-1;
        out(i, j) = dist;
        out(j, i) = dist;
      }else if(size == 9){
        std::vector<int> x = transform_9(std::make_pair(c[i], c[j]));
        long long x_packed = pack_single(x, 1);
        int dist = distmaps::read_pair(x_packed, 1)-1;
        out(i, j) = dist;
        out(j, i) = dist;
      }
    }
  }
  return out;
}*/

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix retrieve_dist_pair(std::vector<std::vector<int> > c1, std::vector<std::vector<int> > c2, int size){ 
  if(size != 13 && size != 9){
    // only support distance pairs for size 9 and 13
    cerr << __FUNCTION__ << ": pairwise distance supported only for cassettes of size 13, 9" << endl;
    Rcpp::stop("Not supported");
  }
  int size_idx = get_size_idx(size);
  Rcpp::NumericMatrix out(c1.size(), c2.size());
  for(int i = 0; i < c1.size(); i++){
    for(int j = 0; j < c2.size(); j++){
      if(size == 13){
        std::vector<int> x = transform_13(std::make_pair(c1[i], c2[j]));
        long long x_packed = pack_single(x, 0);
        int dist = distmaps::read_pair(x_packed, 0)-1;
        out(i, j) = dist;
      }else if(size == 9){
        std::vector<int> x = transform_9(std::make_pair(c1[i], c2[j]));
        long long x_packed = pack_single(x, 1);
        int dist = distmaps::read_pair(x_packed, 1)-1;
        out(i, j) = dist;
      }
    }
  }
  return out;
}

