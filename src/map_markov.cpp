/* implementation of Markov chain model
 * for cassette generation using std::map
 */

#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "pack.h"
#include "cassutil.h"
#include "dist.h"
#include "tables.h"

/* following implementations taken from
 * wehi-project-19/markov 
 **/ 

void gen_odd_recom(const int* o, int* r, int size_idx){
  /** 
   * o is the starting cassette
   * r will contain the list of output ID's
   * size_idx specifies the size of o
   */
  int temp[sizes_odd[size_idx]];
  for(int i = 0; i < rec_counts[size_idx]; i++){
    int* c = odd_rec[size_idx][i];
    for(int j = 0; j < rec_sizes_odd[size_idx][i]; j++){
      // carry out transformation into temp
      temp[j] = o[c[j]];
    }
    int rec_size = rec_sizes_odd[size_idx][i] + rec_sizes_even[size_idx][i];
    // if(rec_size == 1) continue;
    int rec_size_idx = size_to_size_idx2[rec_size];
    if(rec_size_idx == 5){
      r[i] = temp[0];
    }else{
      r[i] = pack_odd(temp, rec_size_idx);
    } 
  }
//  cout << "gen_odd exited happily" << endl;
}

void gen_even_recom(const int* e, int* r, int size_idx){
  /**
   * e is the starting cassette
   * r will contain the list of output ID's
   * s will contain the list of respective size indices
   * size_idx specifies the size of e
   */
  int temp[sizes_even[size_idx]];
  for(int i = 0; i < rec_counts[size_idx]; i++){
    int* c = even_rec[size_idx][i];
    for(int j = 0; j < rec_sizes_even[size_idx][i]; j++){
      // carry out transformation into temp
      temp[j] = e[c[j]];
    }
    int rec_size = rec_sizes_odd[size_idx][i] + rec_sizes_even[size_idx][i];
    // if(rec_size == 1) continue;
    int rec_size_idx = size_to_size_idx2[rec_size];
    if(rec_size_idx == 5){
      // do nothing 
      r[i] = 0;
    }else{
      r[i] = pack_even(temp, rec_size_idx);
    }
  }
//  cout << "gen_even exited happily" << endl;
}

void gen_sign_recom(const int* s, int* r, int size_idx){
    int temp[sizes_odd[size_idx] + sizes_even[size_idx]];
    for(int i = 0; i < rec_counts[size_idx]; i++){
      int* c = sign_rec[size_idx][i];
      for(int j = 0; j < sizes_odd[size_idx] + sizes_even[size_idx]; j++){
        temp[j] = c[j] < 0 ? s[abs(c[j])-1] ^ 1 : s[abs(c[j])-1]; 
      }
      int rec_size = rec_sizes_odd[size_idx][i] + rec_sizes_even[size_idx][i];
      // if(rec_size == 1) continue;
      int rec_size_idx = size_to_size_idx2[rec_size];
      if(rec_size_idx == 5){
        r[i] = temp[0];
      }else{
        r[i] = pack_sign(temp, rec_size_idx);
      }
    }
//    cout << "gen_sign exited happily" << endl;
}

void get_even(int* c, int* e, int size_idx){
    for(int i = 0; i < sizes_even[size_idx]; i++){
        e[i]=( (abs(c[2*i+1]) >> 1)-1 );
    }
}

void get_odd(int* c, int* o, int size_idx){
    for(int i = 0; i < sizes_odd[size_idx]; i++){
        o[i] = ( (abs(c[2*i]) - 1) >> 1 );
    }
}

void get_sign(int* c, int* s, int size_idx){
    /* sign is 1 for reverse */
    for(int i = 0; i < sizes_even[size_idx] + sizes_odd[size_idx]; i++){
        s[i] = ((c[i] < 0) ? 1 : 0);
    }
}

double tr(int size_idx){
  return p[size_idx];
}

//' @export
// [[Rcpp::export]]
Rcpp::List run_markov(int n){
    // sizes 13, 9, 7, 5, 3, 1
    // std::map of [cassette id:prob]
    // outer std::vector corresponds to step
    // inner std::vector corresponds to size_idx
    std::vector<std::vector<std::map<long long, double> > > p(n, std::vector<std::map<long long, double> >(6));
    p[0][0][0] = 1;
    for(int i = 0; i < n-1; i++){
      // i is the previous state
      // i+1 is the next state
      for(int size_idx = 0; size_idx < 5; size_idx++){
        int o[sizes_odd[size_idx]], e[sizes_even[size_idx]], s[sizes_odd[size_idx]+sizes_even[size_idx]];
        int o_rec[rec_counts[size_idx]], e_rec[rec_counts[size_idx]], s_rec[rec_counts[size_idx]];
        // p[i][size_idx] is prev
        // p[i+1][size_idx] is next
        for(auto it = p[i][size_idx].begin(); it != p[i][size_idx].end(); ++it){
          // iterate over cassettes in prev
          unpack(it->first, o, e, s, size_idx);
          gen_odd_recom(o, o_rec, size_idx); gen_even_recom(e, e_rec, size_idx); gen_sign_recom(s, s_rec, size_idx);

          for(int j = 0; j < rec_counts[size_idx]; j++){
            int rec_size_idx = size_to_size_idx2[rec_sizes_odd[size_idx][j] + rec_sizes_even[size_idx][j]]; 
            long long new_cassette;
            if(rec_size_idx == 5){
              // 1-element cassette
              new_cassette = o_rec[j] + 7*s_rec[j];
            }else{
              new_cassette = pack(o_rec[j], e_rec[j], s_rec[j], rec_size_idx);
            }
            p[i+1][rec_size_idx][new_cassette] += it->second*tr(size_idx);
          }
        }
      }
      // special case of size_idx == 5, we have single element cassettes
      int size_idx = 5; 
      for(auto it = p[i][size_idx].begin(); it != p[i][size_idx].end(); ++it){
        p[i+1][size_idx][it->first] += it->second;
      }
    }
    /*
    for(int i = 0; i < n; i++){
      cout << "i = " << i << " " << p[i][0].size() << endl;
      double sum = 0;
      for(int size_idx = 0; size_idx < 6; size_idx++){
        for(auto it = p[i][size_idx].begin(); it != p[i][size_idx].end(); ++it){
          sum += it->second;
        }
      }
      cout << "sum: " << sum << endl;
    } */

    Rcpp::List out = Rcpp::List::create();

    for(int i = 0; i < n; i++){
      std::vector<long long> id;
      std::vector<double> prob;
      std::vector<int> size;
      std::vector<std::string> str;
      for(int size_idx = 0; size_idx < 6; size_idx++){
        for(auto it = p[i][size_idx].begin(); it != p[i][size_idx].end(); ++it){
          id.push_back(it->first);
          prob.push_back(it->second);
          size.push_back(size_idx_to_size[size_idx]);
          if(size_idx != 5){
            std::vector<int> k = unpack_to_vec(it->first, size_idx_to_size[size_idx]);
            str.push_back(get_cass_str(k));
          }else{
            int sign = (it->first)/7;
            int odd = (it->first)%7;
            str.push_back(std::to_string((2*odd+1)*(sign == 1 ? -1 : 1)));
          }
        }
      }
      out.push_back(Rcpp::DataFrame::create(Rcpp::Named("id") = id, Rcpp::Named("p") = prob, Rcpp::Named("size")=size, Rcpp::Named("code")=str));
    }
    return out;
}
