#ifndef __BLINK_MAPPING__H_
#define __BLINK_MAPPING__H_

#include <astroio.hpp>
#include <vector>
#include <string>
#include <unordered_map>


class pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const{
    return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
  }
};


std::vector<int> get_pfb_mapping();

std::unordered_map<std::pair<int, int>, std::pair<int, int>, pairhash> get_visibilities_mapping(const std::string& metafits_filename);

template <typename T>
Visibilities<T> reorder_visibilities(const Visibilities<T>& vis, const std::vector<int>& mapping){



}



#endif
