#ifndef CASSUTIL_H
#define CASSUTIL_H

#include <vector>
#include <string>
#include <Rcpp.h>

std::vector<std::vector<int> > get_cass_vec(std::vector<std::string> c);
std::vector<int> unpack_to_vec(long long c, int size);
std::string get_cass_str(std::vector<int> c);
#endif
