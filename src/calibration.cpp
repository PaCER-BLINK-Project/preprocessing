#include "calibration.hpp"
#include <fstream>



void read_solution(std::string filename, solution_t& sol){
    std::ifstream inp {filename,  std::ios::in | std::ios::binary};
    inp.read(reinterpret_cast<char *>(&sol.header), sizeof(header_t));
    if(!inp){
        std::cerr << "read_solution: error while reading input stream." << std::endl;
        throw std::exception();
    }
    size_t nSolutions {sol.header.intervalCount * sol.header.antennaCount * sol.header.channelCount};
    sol.data = new JonesMatrix<double>[nSolutions];
    inp.read(reinterpret_cast<char *>(sol.data), sizeof(JonesMatrix<double>) * nSolutions);
    if(!inp){
        std::cerr << "read_solution: error while reading input stream." << std::endl;
        throw std::exception();
    }
}