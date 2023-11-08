#include <iostream>
#include <cstring>
#include<tuple>
#include "common.hpp"
#include "../src/mapping.hpp"


std::string data_root_dir;

std::vector<int> get_pfb_mapping();

void test_pfb_mapping(){
    auto mapping = get_pfb_mapping();
    if(mapping[56] != 14) throw TestFailed("'test_pfb_mapping' failed.");
    std::cout << "'test_pfb_mapping' passed." << std::endl;
}


void test_visibilities_mapping(){
    const std::string metadata_file {data_root_dir + "/mwa/1276619416/20200619163000.metafits"}; 
	auto mapping = get_visibilities_mapping(metadata_file);
    int ant, pol;
    const int idx {mapping[2]};
    ant = idx / 2;
    pol = idx % 2;
    if(ant != 66 || pol != 1)
	    throw TestFailed("'test_antenna_mapping' failed.");
    std::cout << "'test_antenna_mapping' passed." << std::endl;
}

void test_reordering(){
	const std::string metadata_file {data_root_dir + "/mwa/1276619416/20200619163000.metafits"}; 
	const std::string vis_file {data_root_dir + "/mwa/1276619416/visibilities/1276619416_20200619163000_gpubox24_00.fits"};
    auto vis = Visibilities::from_fits_file(vis_file);
    auto mapping = get_visibilities_mapping(metadata_file);
	auto reord_vis = reorder_visibilities(vis, mapping);
	
    const float expected_values[] {342.5, -416.25, -346., 840., -308.5, 507., 292., 220.75};
    const int n_baselines = 128 / 2 * 129;
    const int matrix_size = n_baselines * 4 * 2;
    const int baseline_idx = 1;
    const int channel = 1;
    const float* computed_values {
        reinterpret_cast<float*>(reord_vis.data()) + channel * matrix_size + baseline_idx * 8 
    };
    for(int i {0}; i < 8; i++){
        if(std::abs(expected_values[i] - computed_values[i]) > 1e-4){
            std::cerr << "expected_value[i] = " << expected_values[i] << ", computed_value[i] = " << computed_values[i] << std::endl;
            throw TestFailed("'test_reordering': values are different.");
        }
    }
    std::cout << "'test_reordering' passed." << std::endl;
}


int main(void){
    char *pathToData {std::getenv(ENV_DATA_ROOT_DIR)};
    if(!pathToData){
        std::cerr << "'" << ENV_DATA_ROOT_DIR << "' environment variable is not set." << std::endl;
        return -1;
    }
    data_root_dir = std::string {pathToData};
    try{
        test_pfb_mapping();
        test_visibilities_mapping();
        test_reordering();

    } catch (TestFailed ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
