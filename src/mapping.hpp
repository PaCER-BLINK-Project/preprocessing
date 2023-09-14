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

template <typename T>
Visibilities<T> reorder_visibilities(const Visibilities<T>& vis, const antenna_mapping_t& mapping){
    // TODO enforce this function to be applied only to MWA legacy data
    const int n_antennas {128};
    const int n_pols {2};
    const int n_channels {static_cast<int>(vis.nFrequencies)};
    const int n_baselines {(n_antennas + 1) * (n_antennas / 2)};
    const int matrix_size {n_baselines * n_pols * n_pols * 2};
    const int n_intervals {static_cast<int>(vis.integration_intervals())};

    Visibilities<T> final_vis {vis};

    for(int interval {0}; interval < n_intervals; interval++){
        T* new_vis {reinterpret_cast<T*>(final_vis.data) + interval * matrix_size * n_channels};
        T* curr_vis {reinterpret_cast<T*>(vis.data) + interval * matrix_size * n_channels};
        for(int ch {0}; ch < n_channels; ch++){
            for(int src_baseline {0}; src_baseline < n_baselines; src_baseline++){
                int source_a1 {static_cast<int>(-0.5 + std::sqrt(0.25 + 2*src_baseline))};
                int source_a2 {src_baseline - ((source_a1 + 1) * source_a1)/2};
                for(int src_pol1 {0}; src_pol1 < 2; src_pol1++){
                    for(int src_pol2 {0}; src_pol2 < 2; src_pol2++){
                        int dest_a1, dest_p1, dest_a2, dest_p2;

                        std::tie(dest_a1, dest_p1) = mapping.at({source_a1, src_pol1});
                        std::tie(dest_a2, dest_p2) = mapping.at({source_a2, src_pol2});
                        
                        const int source1_idx {source_a1 * 2 + src_pol1};
                        const int source2_idx {source_a2 * 2 + src_pol2};

                        const int dest1_idx {dest_a1 * 2 + dest_p1};
                        const int dest2_idx {dest_a2 * 2 + dest_p2};
                        bool to_conjugate = !((source1_idx < source2_idx && dest1_idx < dest2_idx) || (source1_idx > source2_idx && dest1_idx > dest2_idx));

                        const int src_idx { ch * matrix_size + src_baseline * 8 +  src_pol1 * 4 + src_pol2 * 2};
                        int new_baseline;
                        int new_idx;
                        if(dest_a2 > dest_a1){
                            // we need to transpose polarisations
                            new_baseline =  get_baseline_from_antenna_pair(dest_a2, dest_a1);
                            new_idx = ch * matrix_size + new_baseline * 8 + dest_p2 * 4 + dest_p1 * 2;
                        }else{
                            new_baseline =  get_baseline_from_antenna_pair(dest_a1, dest_a2);
                            new_idx = ch * matrix_size + new_baseline * 8 + dest_p1 * 4 + dest_p2 * 2;
                        }
                        // copy the real part
                        new_vis[new_idx] = curr_vis[src_idx];
                        // now the imaginary part, and correct if necessary
                        new_vis[new_idx + 1] = to_conjugate ? (-curr_vis[src_idx + 1]) : curr_vis[src_idx + 1];
                    }
                }
            }
        }
    }
    return final_vis;
}



#endif
