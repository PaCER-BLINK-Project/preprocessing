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


using antenna_mapping_t = std::unordered_map<std::pair<int, int>, std::pair<int, int>, pairhash>;

std::vector<int> get_pfb_mapping();

antenna_mapping_t get_visibilities_mapping(const std::string& metafits_filename);


inline int get_baseline_from_antenna_pair(int row, int col){
    return (row * (row + 1)) / 2 + col;
}


Visibilities reorder_visibilities(const Visibilities& vis, const antenna_mapping_t& mapping);



#endif
