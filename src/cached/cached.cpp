#include "cached/cached.h"

namespace cache {

void CachedFunction1D::computeAndSave() {
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

void CachedFunction2D::computeAndSave() {
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

void CachedFunction3D::computeAndSave() {
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

}  // namespace cache