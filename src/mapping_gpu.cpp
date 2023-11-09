#include <astroio.hpp>
#include "mapping.hpp"
#include <exception>


__device__ inline int get_baseline_from_antenna_pair(int row, int col){
    return (row * (row + 1)) / 2 + col;
}

__global__ void reorder_visibilities_kernel(const float *current_vis, const int *mapping, float *new_vis){
    
    const unsigned int n_antennas {128};
    const unsigned int n_pols {2};
    const unsigned int n_baselines {(n_antennas + 1) * (n_antennas / 2)};
    const unsigned int matrix_size {n_baselines * n_pols * n_pols * 2};
    const unsigned int n_channels {gridDim.y};
    const unsigned int interval {blockIdx.x};
    const unsigned int ch {blockIdx.y};

    new_vis = new_vis + interval * matrix_size * n_channels;
    current_vis = current_vis + interval * matrix_size * n_channels;

    for(unsigned int src_baseline {threadIdx.x}; src_baseline < n_baselines; src_baseline += blockDim.x){
        const unsigned int source_a1 {static_cast<unsigned int>(-0.5f + std::sqrt(0.25f + 2*src_baseline))};
        const unsigned int source_a2 {src_baseline - ((source_a1 + 1) * source_a1)/2};
        for(unsigned int src_pol1 {0}; src_pol1 < 2; src_pol1++){
            for(unsigned int src_pol2 {0}; src_pol2 < 2; src_pol2++){
                
                const unsigned int source1_idx {source_a1 * 2 + src_pol1};
                const unsigned int source2_idx {source_a2 * 2 + src_pol2};
                const unsigned int dest1_idx {static_cast<unsigned int>(mapping[source1_idx])};
                const unsigned int dest2_idx {static_cast<unsigned int>(mapping[source2_idx])};
                
                unsigned int dest_a1, dest_p1, dest_a2, dest_p2;
                dest_a1 = dest1_idx / 2;
                dest_p1 = dest1_idx % 2;
                dest_a2 = dest2_idx / 2;
                dest_p2 = dest2_idx % 2;
                
                bool to_conjugate = !((source1_idx < source2_idx && dest1_idx < dest2_idx) || (source1_idx > source2_idx && dest1_idx > dest2_idx));

                const unsigned int src_idx { ch * matrix_size + src_baseline * 8 +  src_pol1 * 4 + src_pol2 * 2};
                unsigned int new_baseline;
                unsigned int new_idx;
                if(dest_a2 > dest_a1){
                    // we need to transpose polarisations
                    new_baseline =  get_baseline_from_antenna_pair(dest_a2, dest_a1);
                    new_idx = ch * matrix_size + new_baseline * 8 + dest_p2 * 4 + dest_p1 * 2;
                }else{
                    new_baseline =  get_baseline_from_antenna_pair(dest_a1, dest_a2);
                    new_idx = ch * matrix_size + new_baseline * 8 + dest_p1 * 4 + dest_p2 * 2;
                }
                // copy the real part
                new_vis[new_idx] = current_vis[src_idx];
                // now the imaginary part, and correct if necessary
                new_vis[new_idx + 1] = to_conjugate ? (-current_vis[src_idx + 1]) : current_vis[src_idx + 1];
            }
        }
    }
}



Visibilities reorder_visibilities_gpu(const Visibilities& vis, const MemoryBuffer<int>& mapping){
    // TODO enforce this function to be applied only to MWA legacy data
    if(num_available_gpus() == 0) throw std::runtime_error("reorder_visibilities_gpu: no GPUs detected.");
    if(!vis.on_gpu()) throw std::runtime_error("reorder_visibilities_gpu: 'vis' is not allocated on GPU.");
    if(!mapping.on_gpu()) throw std::runtime_error("reorder_visibilities_gpu: 'mapping' is not allocated on GPU.");
    
    const unsigned int n_antennas {128u};
    const unsigned int n_pols {2u};
    const unsigned int n_channels {vis.nFrequencies};
    const unsigned int n_baselines {(n_antennas + 1) * (n_antennas / 2)};
    const unsigned int matrix_size {n_baselines * n_pols * n_pols * 2};
    const unsigned int n_intervals {static_cast<unsigned int>(vis.integration_intervals())};

    // TODO: avoid memory copying.
    Visibilities final_vis {vis};

    float* new_vis {reinterpret_cast<float*>(final_vis.data())};
    const float* curr_vis {reinterpret_cast<const float*>(vis.data())};

    const int threads_per_block {1024};
    const dim3 n_blocks {n_intervals, n_channels};

    reorder_visibilities_kernel<<<n_blocks, threads_per_block>>>(curr_vis, mapping.data(), new_vis);
    gpuCheckLastError();
    gpuDeviceSynchronize();
    
    return final_vis;
}