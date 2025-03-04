#ifndef GNUPROP_CACHED_H
#define GNUPROP_CACHED_H

#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "simprop.h"

namespace cache {

// Base class for chaced functions
class CachedFunction {
 protected:
  std::string filename;                     // Output file to save cached data
  std::vector<double> results;              // Results of the computation
  std::vector<std::vector<double>> params;  // Parameters of the axes

 public:
  // Constructor with output file
  explicit CachedFunction(const std::string& outputFile) : filename(outputFile) {
    LOGI << "Starting computation cached fuction in " << filename;
  }

  // Utility to save data to a binary file
  template <typename T>
  void saveToBinaryFile(const std::vector<T>& data,
                        const std::vector<std::vector<double>>& params) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
      throw std::runtime_error("Unable to open file: " + filename);
    }

    for (const auto& param : params) {
      file.write(reinterpret_cast<const char*>(&param[0]), sizeof(T));
      file.write(reinterpret_cast<const char*>(&param[1]), sizeof(T));
      file.write(reinterpret_cast<const char*>(&param[2]), sizeof(T));
    }

    for (const auto& value : data) {
      file.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    file.close();
    LOGW << "Data saved in binary format to " << filename;
  }

  // Pure virtual function for tabulation logic
  virtual void computeAndSave() = 0;

  virtual ~CachedFunction() { LOGI << "Finished computation cached fuction in " << filename; };
};

// Derived class for 1D caching
class CachedFunction1D : public CachedFunction {
 private:
  std::function<double(double)> func;
  std::vector<double> xAxis;

 public:
  CachedFunction1D(const std::string& outputFile, std::function<double(double)> f,
                   std::vector<double> x)
      : CachedFunction(outputFile), func(f), xAxis(x) {}

  void computeAndSave() override;
};

// Derived class for 2D caching
class CachedFunction2D : public CachedFunction {
 private:
  std::function<double(double, double)> func;
  std::vector<double> xAxis;
  std::vector<double> yAxis;

 public:
  CachedFunction2D(const std::string& outputFile, std::function<double(double, double)> f,
                   std::vector<double> x, std::vector<double> y)
      : CachedFunction(outputFile), func(f), xAxis(x), yAxis(y) {}

  void computeAndSave() override;
};

// Derived class for 3D caching
class CachedFunction3D : public CachedFunction {
 private:
  std::function<double(double, double, double)> func;
  std::vector<double> xAxis;
  std::vector<double> yAxis;
  std::vector<double> zAxis;

 public:
  CachedFunction3D(const std::string& outputFile, std::function<double(double, double, double)> f,
                   std::vector<double> x, std::vector<double> y, std::vector<double> z)
      : CachedFunction(outputFile), func(f), xAxis(x), yAxis(y), zAxis(z) {}

  void computeAndSave() override;
};

}  // namespace cache

#endif