#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "simprop.h"

// Base class for tabulated functions
class TabulatedFunction {
 protected:
  std::string filename;  // Output file to save tabulated data

 public:
  explicit TabulatedFunction(const std::string& outputFile) : filename(outputFile) {
    LOGI << "Starting computation tabulated fuction in " << filename;
  }

  // Pure virtual function for tabulation logic
  virtual void computeAndSave() = 0;

  // Utility to save data to a binary file (1D data)
  template <typename T>
  void saveToBinaryFile(const std::vector<T>& data) const {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
      throw std::runtime_error("Unable to open file: " + filename);
    }

    for (const auto& value : data) {
      file.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    file.close();
    LOGW << "Data saved in binary format to " << filename;
  }

  virtual ~TabulatedFunction() = default;
};

// Derived class for 1D tabulation
class TabulatedFunction1D : public TabulatedFunction {
 private:
  std::function<double(double)> func;
  std::vector<double> xAxis;

 public:
  TabulatedFunction1D(const std::string& outputFile, std::function<double(double)> f,
                      std::vector<double> x)
      : TabulatedFunction(outputFile), func(f), xAxis(x) {}

  void computeAndSave() override {
    std::vector<double> results;
    auto progressbar = std::make_shared<simprop::utils::ProgressBar>(xAxis.size());
    auto progressbar_mutex = std::make_shared<std::mutex>();
    progressbar->setMutex(progressbar_mutex);
    progressbar->start("Start caching 1D function");
    for (const auto& x : xAxis) {
      progressbar->update();
      results.push_back(func(x));
    }
    saveToBinaryFile(results);
  }
};

// Derived class for 2D tabulation
class TabulatedFunction2D : public TabulatedFunction {
 private:
  std::function<double(double, double)> func;
  std::vector<double> xAxis;
  std::vector<double> yAxis;

 public:
  TabulatedFunction2D(const std::string& outputFile, std::function<double(double, double)> f,
                      std::vector<double> x, std::vector<double> y)
      : TabulatedFunction(outputFile), func(f), xAxis(x), yAxis(y) {}

  void computeAndSave() override {
    std::vector<double> results;
    auto size = xAxis.size() * yAxis.size();
    auto progressbar = std::make_shared<simprop::utils::ProgressBar>(size);
    auto progressbar_mutex = std::make_shared<std::mutex>();
    progressbar->setMutex(progressbar_mutex);
    progressbar->start("Start caching 2D function");
    for (const auto& x : xAxis) {
      for (const auto& y : yAxis) {
        progressbar->update();
        results.push_back(func(x, y));
      }
    }
    saveToBinaryFile(results);
  }
};

// Derived class for 3D tabulation
class TabulatedFunction3D : public TabulatedFunction {
 private:
  std::function<double(double, double, double)> func;
  std::vector<double> xAxis;
  std::vector<double> yAxis;
  std::vector<double> zAxis;

 public:
  TabulatedFunction3D(const std::string& outputFile,
                      std::function<double(double, double, double)> f, std::vector<double> x,
                      std::vector<double> y, std::vector<double> z)
      : TabulatedFunction(outputFile), func(f), xAxis(x), yAxis(y), zAxis(z) {}

  void computeAndSave() override {
    std::vector<double> results;
    auto size = xAxis.size() * yAxis.size() * zAxis.size();
    auto progressbar = std::make_shared<simprop::utils::ProgressBar>(size);
    auto progressbar_mutex = std::make_shared<std::mutex>();
    progressbar->setMutex(progressbar_mutex);
    progressbar->start("Start caching 3D function");
    for (const auto& x : xAxis) {
      for (const auto& y : yAxis) {
        for (const auto& z : zAxis) {
          progressbar->update();
          results.push_back(func(x, y, z));
        }
      }
    }
    saveToBinaryFile(results);
  }
};
