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
#include "dist.h"

using namespace Rcpp;

std::vector<int> get_cass_vec_single(std::string c){
  std::stringstream s(c);
  std::vector<int> out;
  string k = "";
  while(!s.eof()){
    s >> k;
    if(k != "?"){
    	out.push_back(atoi(k.c_str()));
    }else{
      out.push_back(0);
    }
  }
  return out;
}

std::vector<std::string> get_cass_str_vec_single(std::string c){
	std::stringstream s(c);
	std::vector<std::string> out;
	std::string k = "";
	while(!s.eof()){
		s >> k;
		out.push_back(k);
	}
	return out;
}

std::string get_cass_str(std::vector<int> c){
	std::stringstream s;
	for(int i = 0; i < c.size(); i++){
		s << c[i];
    if(i != c.size()-1) s << " ";
	}
	return s.str();
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

//' Convert vector of cassette strings to list of numeric vectors 
//' 
//' @param c vector of cassette strings 
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

//' Check if cassette is valid
//'
//' Checks that cassette is an allowed size, even and odd elements are well-positioned
//' and that there are no repeated elements. 
//' @param c a list of numeric vectors, or a character vector of decoded loxcodes 
//' @return TRUE if the cassette is valid, otherwise FALSE 
//' @export
// [[Rcpp::export]]
std::vector<bool> is_valid(SEXP c){
  switch(TYPEOF(c)){
    case VECSXP:
      return is_valid_impl<std::vector<std::vector<int> > >(as<std::vector<std::vector<int> > >(c));
      break;
    case STRSXP:
      return is_valid_impl<std::vector<std::string> >(as<std::vector<std::string> >(c));
      break;
  }
  return std::vector<bool>();
}

//' Impute missing element in a 13-element cassette
//' @export
// [[Rcpp::export]]
std::string impute_13_impl(std::string c2){
	// impute the missing element
	std::set<int> elems = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
	std::vector<std::string> c = get_cass_str_vec_single(c2);
	std::vector<std::string> out(13);
	std::vector<std::vector<int> > temp(2, std::vector<int>(13));

	if(c.size() != 13){
		// incorrect size 
		return c2; 
	}

	int missing_index = -1;
	for(int i = 0; i < c.size(); i++){
		if(c[i] == "?"){
			// missing index
			missing_index = i;
		}else{
			int k = atoi(c[i].c_str()); // cassette element
			elems.erase(abs(k));
			temp[0][i] = k;
		}
	}
	if(elems.size() != 1 || missing_index < 0){
		// expect only a single unresolved element
		return c2; 
	}
	int missing_elem = *(elems.begin());
	temp[1] = temp[0];
	temp[0][missing_index] = missing_elem;
	temp[1][missing_index] = -missing_elem;

	std::vector<bool> temp_valid = is_valid_impl<std::vector<std::vector<int> > >(temp);
	
	std::vector<long long> packed = pack_impl<std::vector<std::vector<int> > >(temp, temp_valid);

	std::vector<int> dist_orig = retrieve_dist_origin(packed, std::vector<int>({13, 13}));

	//cout << dist_orig[0] << " " << dist_orig[1] << endl;

	int best = dist_orig[0] < dist_orig[1] ? 0 : 1;
	
	return get_cass_str(temp[best]);

}

//' Impute missing elements
//' @export
// [[Rcpp::export]]
std::vector<std::string> impute_13(std::vector<std::string> c, std::vector<int> sizes){
	for(int i = 0; i < c.size(); i++){
		if(sizes[i] == 13){
			c[i] = impute_13_impl(c[i]);
		}
	}
	return c;
}
