#include <astroio.hpp>
#include <metafits_mapping.hpp>
#include "mapping.hpp"
#include <unordered_map>



std::vector<int> get_pfb_mapping(){
    int single_pfb_output_to_input[64] {0, 16, 32, 48};
    for (int i {4}; i < 64; i++){
        single_pfb_output_to_input[i] = single_pfb_output_to_input[i - 4] + 1;
    }
    
    std::vector<int> mapping (256, 0);
    
    for(int p {0}; p < 4; p++){
        for(int i {0}; i < 64; i++){
            mapping[p * 64 + i] = single_pfb_output_to_input[i] + p*64;
        }
    }
    return mapping;
}

std::unordered_map<std::pair<int, int>, std::pair<int, int>, pairhash> get_visibilities_mapping(const std::string& metafits_filename){
    auto metafits_mapping = read_metafits_mapping(metafits_filename);
    auto pfb_mapping = get_pfb_mapping();
    std::unordered_map<std::pair<int, int>, std::pair<int, int>, pairhash> mapping;
    
    for(int source_ant {0}; source_ant < 128; source_ant++){
        for(int source_pol {0}; source_pol < 2; source_pol++){
            int pfb_pos = pfb_mapping[source_ant * 2 + source_pol];
            int final_pos = metafits_mapping[pfb_pos];
            int actual_ant = final_pos / 2;
            int actual_pol = final_pos % 2;
            mapping.insert({{source_ant, source_pol}, {actual_ant, actual_pol}});
        }
    }
    return mapping;
}



Visibilities reorder_visibilities(const Visibilities& vis, const antenna_mapping_t& mapping){
    // TODO enforce this function to be applied only to MWA legacy data
    const int n_antennas {128};
    const int n_pols {2};
    const int n_channels {static_cast<int>(vis.nFrequencies)};
    const int n_baselines {(n_antennas + 1) * (n_antennas / 2)};
    const int matrix_size {n_baselines * n_pols * n_pols * 2};
    const int n_intervals {static_cast<int>(vis.integration_intervals())};

    Visibilities final_vis {vis};

    for(int interval {0}; interval < n_intervals; interval++){
        float* new_vis {reinterpret_cast<float*>(final_vis.data()) + interval * matrix_size * n_channels};
        const float* curr_vis {reinterpret_cast<const float*>(vis.data()) + interval * matrix_size * n_channels};
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