#ifndef __BLINK_MAPPING__H_
#define __BLINK_MAPPING__H_

#include <astroio.hpp>
#include <string>


MemoryBuffer<int> get_visibilities_mapping(const std::string& metafits_filename);


/**
 * @brief Reorder visibilities according to the input-to-antenna index mapping.
 * @param vis: Visibilities object containing visibilities to reorder.
 * @param mapping: MemoryBuffer object wrapping an array of integers, representing
 * the mapping between correlator inputs and antenna indexes.
 * @details For historical reasons the input given to the PFB and then correlator 
 * 
 * @return A Visibilities object containing reordered data.
*/
Visibilities reorder_visibilities(const Visibilities& vis, const MemoryBuffer<int>& mapping);

/**
 * @brief CPU implementation of `reorder_visibilities`.
*/
Visibilities reorder_visibilities_cpu(const Visibilities& vis, const MemoryBuffer<int>& mapping);

#endif
