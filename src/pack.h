#ifndef PACK_H
#define PACK_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <numeric>
using namespace std;

extern const int sizes_odd[5]; 
extern const int sizes_even[5];

void fill_tables();
long long int pack(int o, int e, int s, int size_idx);
void unpack(long long int index, int* o, int* e, int* s, int size_idx);
void join_oes(int* o, int* e, int* s, int* c, int size_idx);
int pack_odd(int* o, int size_idx);
int pack_even(int* e, int size_idx);
int pack_sign(int* s, int size_idx);
unsigned long long xform(unsigned long long a, unsigned long long b);
#endif
