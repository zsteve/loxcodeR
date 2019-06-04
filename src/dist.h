#ifndef DIST_H
#define DIST_H

template<typename T>
std::vector<long long> pack_impl(T c, std::vector<bool> v);

std::vector<int> retrieve_dist_origin(std::vector<long long> c, std::vector<int> sizes);

#endif 
