#ifndef __BLINK_MAPPING__H_
#define __BLINK_MAPPING__H_

#include <astroio.hpp>
#include <string>


MemoryBuffer<int> get_visibilities_mapping(const std::string& metafits_filename);


Visibilities reorder_visibilities(const Visibilities& vis, const MemoryBuffer<int>& mapping);

#endif
