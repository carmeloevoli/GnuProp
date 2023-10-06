#ifndef BENIAMINO_COSMOLOGY_H
#define BENIAMINO_COSMOLOGY_H

#include "units.h"

namespace beniamino {

inline double E(double z) { return std::sqrt(SI::OmegaL + SI::OmegaM * pow3(1. + z)); }
inline double hubbleRate(double z) { return SI::H0 * E(z); }
inline double dtdz(double z) { return 1. / hubbleRate(z) / (1. + z); }

}  // namespace beniamino

#endif