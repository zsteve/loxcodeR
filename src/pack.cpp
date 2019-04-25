#include <iostream>
#include <algorithm>
#include <vector>
#include <numeric>

#include "pack.h"
  
using namespace std;

//template <typename Iter>
/* boonext(Iter begin, Iter end)
{
  if (begin == end) {return false;} --end;
  if ((*end & 1) == 0){++*end; return true;}
  else{--*end;return next(begin, end);}
} */

// number of permutations possible
// in order: 13, 9, 7, 5, 3 element cassettes

const int nodd[]={5040, 2520, 840, 210, 42}; 
const int neven[]={720, 360, 120, 30, 6}; 
const int nsign[]={8192, 512, 128, 32, 8};

// powers table
const int pow7[]={1,7,49,343,2401,16807,117649};
const int pow6[]={1,6,36,216,1296,7776};
const int pow2[]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096};

// for the below we will have 'ragged' arrays but for simplicity we just allocate
// the maximal value
int table_odd[5][823543]={0}; 
int table_even[5][46656]={0};
int table_odd_reverse[5][5040]={0};
int table_even_reverse[5][720]={0};

// sizes of odd/even parts ettes
const int sizes_odd[5] = {7, 5, 4, 3, 2};
const int sizes_even[5] = {6, 4, 3, 2, 1};

// debug

extern void print_vec(vector<int> v);

unsigned long long xform(unsigned long long a, unsigned long long b){
  // this definition ONLY for 13-element cassesttes (size_idx = 0)
  // take a as the origin and transform b
  int o_a[7], e_a[6], s_a[13];
  int o_b[7], e_b[6], s_b[13];
  unpack(a, o_a, e_a, s_a, 0);
  unpack(b, o_b, e_b, s_b, 0);
  int o_c[7], e_c[6], s_c[13];  // transformed b
  // position indices
  int p_o[7], p_e[6];
  for(int i = 0; i < 7; i++){ // odd position indices
    p_o[o_a[i]] = i;
  }

  for(int i = 0; i < 6; i++){ // even position indices
    p_e[e_a[i]] = i;
  }

  // transform odd element indices and signs
  for(int i = 0; i < 7; i++){
    o_c[i] = p_o[o_b[i]];
    s_c[2*i] = (s_a[2*p_o[o_b[i]]] == s_b[2*i]) ? 0 : 1;
  }


  // transform even element indices and signs
  for(int i = 0; i < 6; i++){
    e_c[i] = p_e[e_b[i]];
    s_c[2*i + 1] = (s_a[2*p_e[e_b[i]]+1] == s_b[2*i + 1]) ? 0 : 1;
  }


  unsigned long long t_index = pack(pack_odd(o_c, 0), pack_even(e_c, 0), pack_sign(s_c, 0), 0);

  return t_index;

}

void fill_tables(){
  for(int size_idx = 0; size_idx < 5; ++size_idx){
    vector<int> O_mask(7); std::fill(O_mask.end() - sizes_odd[size_idx], O_mask.end(), 1);
    vector<int> E_mask(6); std::fill(E_mask.end() - sizes_even[size_idx], E_mask.end(), 1); 

    int count=0;
    do{
      vector<int> O(sizes_odd[size_idx]); for(int i = 0, j = 0; i < O_mask.size(); ++i) if(O_mask[i]) O[j++] = i; 
      do{
        table_odd[size_idx][(std::inner_product(std::begin(O), std::end(O), std::begin(pow7), 0))]=count;
        table_odd_reverse[size_idx][count]=(std::inner_product(std::begin(O), std::end(O), std::begin(pow7), 0));
        count++;
      }while ( std::next_permutation(O.begin(),O.end()) );
    }while ( std::next_permutation(O_mask.begin(), O_mask.end()));
  
    count=0;
    do{ 
      vector<int> E(sizes_even[size_idx]); for(int i = 0, j = 0; j < E_mask.size(); ++i) if(E_mask[i]) E[j++] = i;
      do{
        table_even[size_idx][(std::inner_product(std::begin(E), std::end(E), std::begin(pow6), 0))]=count;
        table_even_reverse[size_idx][count]=(std::inner_product(std::begin(E), std::end(E), std::begin(pow6), 0));
        count++;
      }while ( std::next_permutation(E.begin(),E.end()) );
    }while (std::next_permutation(E_mask.begin(), E_mask.end()));
  }
}

int inner_product(const int* a,const int* b, int N){
  int r = 0;
  for(int i = 0; i < N; i++){
    r+=a[i]*b[i];
  }
  return r;
}

long long int pack(int o, int e, int s, int size_idx)
{
  return ((long long)o * neven[size_idx] * nsign[size_idx] + (long long)e * nsign[size_idx] + (long long)s);
};


void join_oes(int* o, int* e, int* s, int* c, int size_idx){
	// join odd-even-sign to form a single cassette
	// size_idx/2 + 1 odd elements
	// size_idx/2 even elements
	for(int i = 0; i < sizes_odd[size_idx] + sizes_even[size_idx]; ++i){
		if(i % 2 == 0){
			// odd
			c[i] = 2*o[i/2]+1;
		}else{
			// even
			c[i] = (e[i/2]+1)*2;
		}
		c[i]*=(s[i] == 1 ? -1 : 1);
	}
	return;
}

int pack_odd(int* o, int size_idx){
  return table_odd[size_idx][inner_product(o, pow7, sizes_odd[size_idx])];
}

int pack_even(int* e, int size_idx){
  return table_even[size_idx][inner_product(e, pow6, sizes_even[size_idx])];
}

int pack_sign(int* s, int size_idx){
  return inner_product(s, pow2, sizes_odd[size_idx] + sizes_even[size_idx]);
}

void unpack(long long int index, int* o, int* e, int* s, int size_idx)
{
  int i = index / (nsign[size_idx] * neven[size_idx]) ;             // o 
  int j = (index % (nsign[size_idx] * neven[size_idx])) / nsign[size_idx];    // e
  int k = index - j * nsign[size_idx] - i * nsign[size_idx] * neven[size_idx];// s 

  i=table_odd_reverse[size_idx][i];
  j=table_even_reverse[size_idx][j];

  for(int a=0;a<sizes_odd[size_idx];a++){div_t  d = div(i,7); o[a]=d.rem; i =d.quot ;}
  for(int a=0;a<sizes_even[size_idx];a++){div_t  d = div(j,6); e[a]=d.rem; j =d.quot ;}
  for(int a=0;a<(sizes_odd[size_idx] + sizes_even[size_idx]);a++){div_t d = div(k,2); s[a]=d.rem; k =d.quot ;}
};
