#include <iostream>
#include <cstring>
#include<tuple>
#include "common.hpp"
#include "../src/mapping.hpp"


std::string data_root_dir;


void test_pfb_mapping(){
    auto mapping = get_pfb_mapping();
    if(mapping[56] != 14) throw TestFailed("'test_pfb_mapping' failed.");
    std::cout << "'test_pfb_mapping' passed." << std::endl;
}


void test_visibilities_mapping(){
    const std::string metadata_file {data_root_dir + "/mwa/1103645160/20200619163000.metafits"}; 
	auto mapping = get_visibilities_mapping(metadata_file);
    int ant, pol;
    std::tie(ant, pol) = mapping[{1, 0}];
    if(ant != 66 || pol != 1)
	    throw TestFailed("'test_antenna_mapping' failed.");
    std::cout << "'test_antenna_mapping' passed." << std::endl;
}

void test_reordering(){
	const std::string metadata_file {data_root_dir + "/mwa/1103645160/20200619163000.metafits"}; 
	
	//auto reord_vis = reorder_visibilities(vis, meta);
	
	throw TestFailed("'test_reordering': Implement me.");
    std::cout << "'test_reordering' passed." << std::endl;
}


int main(void){
    char *pathToData {std::getenv(ENV_DATA_ROOT_DIR)};
    if(!pathToData){
        std::cerr << "'" << ENV_DATA_ROOT_DIR << "' environment variable is not set." << std::endl;
        return -1;
    }
    data_root_dir = {pathToData};
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
