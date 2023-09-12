#include <astroio.hpp>
#include <metafits_mapping.hpp>
#include "mapping.hpp"
#include <unordered_map>

template <>
Visibilities<float> reorder_visibilities(const Visibilities<float>& vis, const std::vector<int>& mapping);


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