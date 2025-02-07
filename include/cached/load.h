#ifndef GNUPROP_CACHE_LOAD_H
#define GNUPROP_CACHE_LOAD_H

#include <string>
#include <vector>

namespace cache {

std::vector<double> loadFromFile1D(const std::string& filename);
std::vector<std::vector<double>> loadFromFile(const std::string& filename, size_t redshiftSize);

}  // namespace cache

#endif