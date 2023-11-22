#ifndef __CALIBRATION_GPU_H__
#define __CALIBRATION_GPU_H__

#include "calibration.hpp"
#include <gpu_macros.hpp>



void apply_solutions_gpu(Visibilities &vis, const CalibrationSolutions& sol, unsigned int coarse_channel_index);

#endif
