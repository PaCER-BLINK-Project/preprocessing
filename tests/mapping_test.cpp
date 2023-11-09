#include <iostream>
#include <cstring>
#include <tuple>
#include <astroio.hpp>
#include "common.hpp"
#include "../src/mapping.hpp"
#ifdef __GPU__
#include "../src/mapping_gpu.hpp"
#endif

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
	auto reord_vis = reorder_visibilities_cpu(vis, mapping);
	
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

#ifdef __GPU__
void test_reordering_gpu(){

    const std::string metadata_file {data_root_dir + "/mwa/1276619416/20200619163000.metafits"}; 
	const std::string vis_file {data_root_dir + "/mwa/1276619416/visibilities/1276619416_20200619163000_gpubox24_00.fits"};
    auto vis = Visibilities::from_fits_file(vis_file);
    auto mapping = get_visibilities_mapping(metadata_file);
	auto reord_vis = reorder_visibilities_cpu(vis, mapping);
    vis.to_gpu();
    mapping.to_gpu();
    auto reord_vis_gpu = reorder_visibilities_gpu(vis, mapping);
    reord_vis_gpu.to_cpu();
    
    if(reord_vis.size() != reord_vis_gpu.size()) throw TestFailed("test_reordering_gpu: lengths differ!");
    for(unsigned long long i {0}; i < reord_vis_gpu.size(); i++){
        if (reord_vis[i].real() != reord_vis_gpu[i].real() ||  reord_vis[i].imag() != reord_vis_gpu[i].imag()){
            throw TestFailed("test_reordering_gpu: elements differ!");
        }
    }
    std::cout << "'test_reordering_gpu' passed." << std::endl;
}
#endif


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
        #ifdef __GPU__
        test_reordering_gpu();
        #endif
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
