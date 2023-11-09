#ifndef __BLINK_MAPPING_GPU__H_
#define __BLINK_MAPPING_GPU__H_

#include <astroio.hpp>


Visibilities reorder_visibilities_gpu(const Visibilities& vis, const MemoryBuffer<int>& mapping);

#endif
