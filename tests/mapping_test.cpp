#include <iostream>
#include <cstring>
#include "common.hpp"
#include "../src/mapping.hpp"


std::string data_root_dir;


void test_reordering(){
	CObsMetadata meta;
	std::string metadata_file {data_root_dir + "/mwa/1103645160/1103645160.metafits"}; 
	meta.ReadMetaData(metadata_file.c_str());

	auto reord_vis = reorder_visibilities(vis, meta);
	
	throw TestFailed("Implement me.");
    std::cout << "'test_reord' passed." << std::endl;
}


int main(void){
    char *pathToData {std::getenv(ENV_DATA_ROOT_DIR)};
    if(!pathToData){
        std::cerr << "'" << ENV_DATA_ROOT_DIR << "' environment variable is not set." << std::endl;
        return -1;
    }
    data_root_dir = {pathToData};
    try{
        
        test_reordering();

    } catch (TestFailed ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
