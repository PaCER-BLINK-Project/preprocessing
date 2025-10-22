#include "calibration.hpp"
#include <memory_buffer.hpp>
#include <fstream>
#ifdef __GPU__
#include "calibration_gpu.hpp"
#endif


CalibrationSolutions CalibrationSolutions::from_file(const std::string& filename){
    std::ifstream inp {filename,  std::ios::in | std::ios::binary};
    CalibrationSolutions::Header header;
    inp.read(reinterpret_cast<char *>(&header), sizeof(CalibrationSolutions::Header));
    if(!inp){
        std::cerr << "read_solution: error while reading input stream." << std::endl;
        throw std::exception();
    }
    size_t n_solutions {header.interval_count * header.antenna_count * header.channel_count};
    MemoryBuffer<JonesMatrix<double>> mb_sol {n_solutions};
    inp.read(reinterpret_cast<char*>(mb_sol.data()), sizeof(JonesMatrix<double>) * n_solutions);
    if(!inp){
        std::cerr << "read_solution: error while reading input stream." << std::endl;
        throw std::exception();
    }
    return {std::move(mb_sol), header};
}



void apply_solutions(Visibilities &vis, CalibrationSolutions& sol, unsigned int coarse_channel_index){
    #ifdef __GPU__
    if(gpu_support() && num_available_gpus() > 0 && vis.on_gpu()){
        sol.to_gpu();
        return apply_solutions_gpu(vis, sol, coarse_channel_index);
    }else{
        sol.to_cpu();
        return apply_solutions_cpu(vis, sol, coarse_channel_index);
    }
    #else
    return apply_solutions_cpu(vis, sol, coarse_channel_index);
    #endif
}



void apply_solutions_cpu(Visibilities &vis, const CalibrationSolutions& sol, unsigned int coarse_channel_index){
    if(sol.header.channel_count % 24 != 0) throw std::exception();
    if(sol.header.antenna_count != vis.obsInfo.nAntennas) throw std::exception();
    
    const unsigned int n_solfreq_per_band = sol.header.channel_count / 24u;
    const unsigned int solutions_offset = coarse_channel_index * n_solfreq_per_band;
    
    int channelRatio;
    if(n_solfreq_per_band > vis.nFrequencies)
        channelRatio = n_solfreq_per_band / vis.nFrequencies;
    else
        channelRatio =  vis.nFrequencies / n_solfreq_per_band;

    const unsigned int n_baselines {((vis.obsInfo.nAntennas + 1) * vis.obsInfo.nAntennas) / 2};
    const size_t matrix_size {n_baselines * vis.obsInfo.nPolarizations * vis.obsInfo.nPolarizations * 2};
    const size_t n_interval_values {matrix_size * vis.nFrequencies};

    #pragma omp parallel for schedule(static) if(vis.integration_intervals() <= 10)
    for(unsigned int nInterval = 0; nInterval < vis.integration_intervals(); nInterval++){
       #pragma omp parallel for schedule(static) if(vis.integration_intervals() <= 10)
        for(unsigned int baseline = 0; baseline < n_baselines; baseline++){
            for(unsigned int ch {0}; ch < vis.nFrequencies; ch++){
                unsigned int solChannel;
                if(n_solfreq_per_band > vis.nFrequencies)
                    solChannel = solutions_offset + ch * channelRatio;
                else
                    solChannel = solutions_offset + ch / channelRatio;
            
                unsigned int a1 {static_cast<unsigned int>(-0.5 + std::sqrt(0.25 + 2*baseline))};
                unsigned int a2 {baseline - ((a1 + 1) * a1)/2};
                const JonesMatrix<double> &solA = sol.data()[a1 * sol.header.channel_count + solChannel];
                const JonesMatrix<double> &solB = sol.data()[a2 * sol.header.channel_count + solChannel];
                if(solA.isnan() || solB.isnan()) continue;
                float *data {reinterpret_cast<float*>(vis.data()) + nInterval * n_interval_values + matrix_size * ch + 8 * baseline};
                
                JonesMatrix<double> visData {JonesMatrix<double>::from_array<double, float>(data)};
                JonesMatrix<double> tmp = calibrate_visibility(visData, solA, solB);

                data[0] = static_cast<float>(tmp.XX.real);
                data[1] = static_cast<float>(tmp.XX.imag);
                data[2] = static_cast<float>(tmp.XY.real);
                data[3] = static_cast<float>(tmp.XY.imag);
                data[4] = static_cast<float>(tmp.YX.real);
                data[5] = static_cast<float>(tmp.YX.imag);
                data[6] = static_cast<float>(tmp.YY.real);
                data[7] = static_cast<float>(tmp.YY.imag);
            }
        }
    }
}
