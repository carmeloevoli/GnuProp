#include "cached/load.h"

#include <cassert>
#include <fstream>

namespace cache {

std::vector<double> loadFromFile1D(const std::string& filename) {
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file: " + filename);
  }

  // Get file size
  file.seekg(0, std::ios::end);
  std::streamsize fileSize = file.tellg();
  file.seekg(0, std::ios::beg);

  // Check if file size is a multiple of double size
  if (fileSize % sizeof(double) != 0) {
    throw std::runtime_error("Invalid binary file size, not a multiple of double.");
  }

  // Calculate the number of doubles
  size_t numDoubles = fileSize / sizeof(double);

  // Read the data into a vector
  std::vector<double> data(numDoubles);
  if (!file.read(reinterpret_cast<char*>(data.data()), fileSize)) {
    throw std::runtime_error("Error reading binary file.");
  }

  file.close();

  return data;
}

std::vector<std::vector<double>> loadFromFile(const std::string& filename, size_t redshiftSize) {
  // Open file
  std::ifstream file(filename, std::ios::binary);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file: " + filename);
  }

  // Get file size
  file.seekg(0, std::ios::end);
  std::streamsize fileSize = file.tellg();
  file.seekg(0, std::ios::beg);

  // Check if file size is a multiple of double size
  if (fileSize % sizeof(double) != 0) {
    throw std::runtime_error("Invalid binary file size, not a multiple of double.");
  }

  // Calculate the number of doubles
  size_t numDoubles = fileSize / sizeof(double);

  // Read the data into a vector
  std::vector<double> data(numDoubles);
  if (!file.read(reinterpret_cast<char*>(data.data()), fileSize)) {
    throw std::runtime_error("Error reading binary file.");
  }

  // Close file
  file.close();

  // Slice data in redshift chunks
  std::vector<std::vector<double>> result;
  size_t sliceSize = numDoubles / redshiftSize;
  assert(sliceSize * redshiftSize == numDoubles);
  for (size_t i = 0; i < redshiftSize; i++) {
    result.emplace_back(data.begin() + (i * sliceSize), data.begin() + ((i + 1) * sliceSize - 1));
  }

  return result;
}

}  // namespace cache