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

// void CachedFunction3D::computeAndSave() {
//   std::vector<double> results;
//   std::mutex results_mutex;  // Mutex for thread-safe access to results

//   auto size = xAxis.size() * yAxis.size() * zAxis.size();
//   auto progressbar = std::make_shared<simprop::utils::ProgressBar>(size);
//   auto progressbar_mutex = std::make_shared<std::mutex>();
//   progressbar->setMutex(progressbar_mutex);
//   progressbar->start("Start caching 3D function");

//   std::vector<std::thread> threads;
//   std::mutex progress_mutex;

//   // Function for each thread
//   auto worker = [&](size_t start, size_t end) {
//     std::vector<double> local_results;  // Thread-local storage to minimize lock contention

//     for (size_t i = start; i < end; ++i) {
//       for (const auto& y : yAxis) {
//         for (const auto& z : zAxis) {
//           local_results.push_back(func(xAxis[i], y, z));

//           // Update progress bar safely
//           std::lock_guard<std::mutex> lock(progress_mutex);
//           progressbar->update();
//         }
//       }
//     }

//     // Merge results into the main vector safely
//     std::lock_guard<std::mutex> lock(results_mutex);
//     results.insert(results.end(), local_results.begin(), local_results.end());
//   };

//   size_t num_threads = std::thread::hardware_concurrency();
//   size_t chunk_size = xAxis.size() / num_threads;
//   size_t remainder = xAxis.size() % num_threads;

//   size_t start = 0;
//   for (size_t i = 0; i < num_threads; ++i) {
//     size_t end = start + chunk_size + (i < remainder ? 1 : 0);
//     if (start < end) {
//       threads.emplace_back(worker, start, end);
//     }
//     start = end;
//   }

//   // Join threads
//   for (auto& t : threads) {
//     t.join();
//   }

//   saveToBinaryFile(results);
// }

// #include <mutex>
// #include <thread>
// #include <vector>

// void CachedFunction3D::computeAndSave() {
//   std::vector<double> results;
//   std::mutex results_mutex;  // Mutex for thread-safe access to results

//   auto size = xAxis.size() * yAxis.size() * zAxis.size();
//   auto progressbar = std::make_shared<simprop::utils::ProgressBar>(size);
//   auto progressbar_mutex = std::make_shared<std::mutex>();
//   progressbar->setMutex(progressbar_mutex);
//   progressbar->start("Start caching 3D function");

//   std::vector<std::thread> threads;
//   std::mutex progress_mutex;

//   // Function for each thread
//   auto worker = [&](size_t start, size_t end) {
//     std::vector<double> local_results;  // Thread-local storage to minimize lock contention

//     for (size_t i = start; i < end; ++i) {
//       for (const auto& y : yAxis) {
//         for (const auto& z : zAxis) {
//           local_results.push_back(func(xAxis[i], y, z));

//           // Update progress bar safely
//           std::lock_guard<std::mutex> lock(progress_mutex);
//           progressbar->update();
//         }
//       }
//     }

//     // Merge results into the main vector safely
//     std::lock_guard<std::mutex> lock(results_mutex);
//     results.insert(results.end(), local_results.begin(), local_results.end());
//   };

//   size_t num_threads = std::thread::hardware_concurrency();
//   size_t chunk_size = xAxis.size() / num_threads;
//   size_t remainder = xAxis.size() % num_threads;

//   size_t start = 0;
//   for (size_t i = 0; i < num_threads; ++i) {
//     size_t end = start + chunk_size + (i < remainder ? 1 : 0);
//     if (start < end) {
//       threads.emplace_back(worker, start, end);
//     }
//     start = end;
//   }

//   // Join threads
//   for (auto& t : threads) {
//     t.join();
//   }

//   saveToBinaryFile(results);
// }

}  // namespace cache