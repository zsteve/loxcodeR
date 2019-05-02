#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <iostream>

using namespace std;

#define SIGN(x) (x > 0 ? 1 : -1)

void  transform_9(std::pair<std::vector<int>, std::vector<int> > c){
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
//  cout << "fwd map" << endl;
//  for(int i = 0; i < 13; i++){
//    cout << elt_map[i+1] << " ";
//  }cout << endl;

  // now flip map
  std::map<int, int> rev_map;
  for(int i = 0; i < 13; i++){
    rev_map[abs(elt_map[i+1])] = SIGN(elt_map[i+1])*(i+1);
  }
  
//  cout << "rev map" << endl;
//  for(int i = 0; i < 13; i++){
//    cout << rev_map[i+1] << " ";
//  }cout << endl;
  vector<int> out(9);
  for(int i = 0; i < 9; i++){
    out[i] = rev_map[abs(c.second[i])]*SIGN(c.second[i]);
  }
  for(int i = 0; i < 9; i++) cout << out[i] << " "; cout << endl;
  return;
}
int main(){
  // std::vector<int> a = {5, 6, 7, 8, 9, 10, 11, 12, 13};
  std::vector<int> a = {1, -2, 5, 8, 7, -4, 3, 10, -9};
  std::vector<int> b = {5, -6, 7, -8, 9, -10, 11, -12, 13} ;
  transform_9(std::make_pair(a, b));

}
