#include "calibration.hpp"
#include <memory_buffer.hpp>
#include <fstream>

__global__ void apply_solutions_kernel(float* vis, ObservationInfo obsInfo, const JonesMatrix<double>* const sol,
        CalibrationSolutions::Header sol_header, unsigned int coarse_channel_index){

    const unsigned int n_solfreq_per_band = sol_header.channel_count / 24u;
    const unsigned int solutions_offset = coarse_channel_index * n_solfreq_per_band;
    const unsigned int n_baselines {(obsInfo.nAntennas + 1) * (obsInfo.nAntennas / 2)};
    const size_t matrix_size {n_baselines * obsInfo.nPolarizations * obsInfo.nPolarizations * 2};
    const unsigned int n_channels {gridDim.y};
    const size_t n_interval_values {matrix_size * n_channels};
    const unsigned int interval {blockIdx.x};
    const unsigned int ch {blockIdx.y};

    int channelRatio;
    if(n_solfreq_per_band > n_channels)
        channelRatio = n_solfreq_per_band / n_channels;
    else
        channelRatio =  n_channels / n_solfreq_per_band;

    for(unsigned int baseline {threadIdx.x}; baseline < n_baselines; baseline+= blockDim.x){
        unsigned int solChannel;
        if(n_solfreq_per_band > n_channels)
            solChannel = (ch + solutions_offset) * channelRatio;
        else
            solChannel = (ch + solutions_offset) / channelRatio;
    
        unsigned int a1 {static_cast<unsigned int>(-0.5 + std::sqrt(0.25 + 2*baseline))};
        unsigned int a2 {baseline - ((a1 + 1) * a1)/2};
        const JonesMatrix<double> &solA = sol[a1 * sol_header.channel_count + solChannel];
        const JonesMatrix<double> &solB = sol[a2 * sol_header.channel_count + solChannel];
        float *data {vis + interval * n_interval_values + matrix_size * ch + 8 * baseline};
        
        JonesMatrix<double> visData {JonesMatrix<double>::from_array<double, float>(data)};
        JonesMatrix<double> tmp = solA * (visData * solB.conjtrans());

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


void apply_solutions_gpu(Visibilities &vis, const CalibrationSolutions& sol, unsigned int coarse_channel_index){    
    if(sol.header.channel_count % 24 != 0) throw std::exception();
    if(sol.header.antenna_count != vis.obsInfo.nAntennas) throw std::exception();
    if(!vis.on_gpu())  throw std::runtime_error("apply_solutions_gpu: 'vis' is not allocated on GPU.");
    if(!sol.on_gpu())  throw std::runtime_error("apply_solutions_gpu: 'sol' is not allocated on GPU.");

    const int threads_per_block {1024};
    const dim3 n_blocks {static_cast<unsigned int>(vis.integration_intervals()), vis.nFrequencies};

    apply_solutions_kernel<<<n_blocks, threads_per_block>>>(
        reinterpret_cast<float*>(vis.data()), vis.obsInfo, sol.data(), sol.header, coarse_channel_index);
    gpuCheckLastError();
    gpuDeviceSynchronize();
}
