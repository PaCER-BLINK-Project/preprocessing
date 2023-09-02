#include <astroio.hpp>
#include "mapping.hpp"
template <>
Visibilities<float> reorder_visibilities(const Visibilities<float>& vis, CObsMetadata);

