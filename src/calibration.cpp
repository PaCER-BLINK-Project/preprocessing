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



void apply_solution(Visibilities &vis, solution_t& sol, unsigned int coarse_channel_index){
    if(sol.header.channelCount % 24 != 0) throw std::exception();
    if(sol.header.antennaCount != vis.obsInfo.nAntennas) throw std::exception();
    
    const unsigned int n_solfreq_per_band = sol.header.channelCount / 24u;
    const unsigned int solutions_offset = coarse_channel_index * n_solfreq_per_band;
    
    int channelRatio;
    if(n_solfreq_per_band > vis.nFrequencies)
        channelRatio = n_solfreq_per_band / vis.nFrequencies;
    else
        channelRatio =  vis.nFrequencies / n_solfreq_per_band;

    const unsigned int n_baselines {(vis.obsInfo.nAntennas + 1) * (vis.obsInfo.nAntennas / 2)};
    const size_t matrix_size {n_baselines * vis.obsInfo.nPolarizations * vis.obsInfo.nPolarizations * 2};
    const size_t n_interval_values {matrix_size * vis.nFrequencies};

    #pragma omp parallel for schedule(static) if(vis.integration_intervals() <= 10)
    for(unsigned int nInterval = 0; nInterval < vis.integration_intervals(); nInterval++){
       #pragma omp parallel for schedule(static) if(vis.integration_intervals() <= 10)
        for(unsigned int baseline = 0; baseline < n_baselines; baseline++){
            for(unsigned int ch {0}; ch < vis.nFrequencies; ch++){
                unsigned int solChannel;
                if(n_solfreq_per_band > vis.nFrequencies)
                    solChannel = (ch + solutions_offset) * channelRatio;
                else
                    solChannel = (ch + solutions_offset) / channelRatio;
            
                unsigned int a1 {static_cast<unsigned int>(-0.5 + std::sqrt(0.25 + 2*baseline))};
                unsigned int a2 {baseline - ((a1 + 1) * a1)/2};
                JonesMatrix<double> &solA = reinterpret_cast<JonesMatrix<double>&>(sol.data[a1 * sol.header.channelCount + solChannel]);
                JonesMatrix<double> &solB = reinterpret_cast<JonesMatrix<double>&>(sol.data[a2 * sol.header.channelCount + solChannel]);
                float *data {reinterpret_cast<float*>(vis.data()) + nInterval * n_interval_values + matrix_size * ch + 8 * baseline};
                
                JonesMatrix<double> visData {JonesMatrix<double>::from_array<double, float>(data)};
                JonesMatrix<double> tmp = calibrate_visibility(visData, solA, solB);

                data[0] = static_cast<float>(tmp.XX.real());
                data[1] = static_cast<float>(tmp.XX.imag());
                data[2] = static_cast<float>(tmp.XY.real());
                data[3] = static_cast<float>(tmp.XY.imag());
                data[4] = static_cast<float>(tmp.YX.real());
                data[5] = static_cast<float>(tmp.YX.imag());
                data[6] = static_cast<float>(tmp.YY.real());
                data[7] = static_cast<float>(tmp.YY.imag());
            }
        }
    }
}
