#ifndef GNUPROP_INTERPOLATE_H
#define GNUPROP_INTERPOLATE_H

#include <cassert>
#include <cmath>
#include <vector>

namespace utils {

inline std::size_t getIndex(double x, double xMin, double xMax, std::size_t N) {
#ifdef DEBUG
  assert(N > 2);
  assert(xMin < xMax);
  assert(x >= xMin && x <= xMax);
#endif

  // Precompute the reciprocal once
  const double invDx = (N - 1) / (xMax - xMin);

  // Floating index in [0, N-1].
  const double floatIndex = (x - xMin) * invDx;
  std::size_t i = static_cast<std::size_t>(std::floor(floatIndex));

  return (i >= N - 1) ? N - 1 : i;
}

inline double interpolate1D(double x, double xMin, double xMax, const std::vector<double>& Y) {
#ifdef DEBUG
  assert(!Y.empty());
  assert(xMin < xMax);
  assert(x >= xMin && x <= xMax);
#endif

  const std::size_t N = Y.size();

  // Precompute the reciprocal once
  const double invDx = (N - 1) / (xMax - xMin);

  // Floating index in [0, N-1].
  const double floatIndex = (x - xMin) * invDx;
  std::size_t i = static_cast<std::size_t>(std::floor(floatIndex));

  // Clamp i so that i+1 doesn't go out of range:
  if (i >= N - 1) {
    // If x is exactly xMax, i can be N-1 => just return last element
    return Y[N - 1];
  }

  // Fractional part
  const double t = floatIndex - static_cast<double>(i);

  // Linear interpolation
  return (1.0 - t) * Y[i] + t * Y[i + 1];
}

// Row-major storage assumed: Z[i, j] = Z[i*ySize + j]
inline double interpolate2D(double x, double y, double xMin, double xMax, std::size_t xSize,
                            double yMin, double yMax, std::size_t ySize,
                            const std::vector<double>& Z) {
#ifdef DEBUG
  assert(xSize > 1 && ySize > 1);
  assert(Z.size() == xSize * ySize);
  assert(xMin < xMax && yMin < yMax);
  assert(x >= xMin && x <= xMax);
  assert(y >= yMin && y <= yMax);
#endif

  // Precompute reciprocals for efficiency
  const double invDx = (xSize - 1) / (xMax - xMin);
  const double invDy = (ySize - 1) / (yMax - yMin);

  // Continuous indices
  const double floatIndexX = (x - xMin) * invDx;
  const double floatIndexY = (y - yMin) * invDy;

  // Integer parts
  std::size_t i = static_cast<std::size_t>(std::floor(floatIndexX));
  std::size_t j = static_cast<std::size_t>(std::floor(floatIndexY));

  // Clamp
  if (i >= xSize - 1) i = xSize - 2;  // if x==xMax => i+1 out of range
  if (j >= ySize - 1) j = ySize - 2;  // likewise for y==yMax

  // Fractional parts
  const double tx = floatIndexX - static_cast<double>(i);
  const double ty = floatIndexY - static_cast<double>(j);

  // Retrieve the four corners
  // Row-major => element (i, j) is at index i*ySize + j
  const double Q11 = Z[i * ySize + j];
  const double Q21 = Z[(i + 1) * ySize + j];
  const double Q12 = Z[i * ySize + (j + 1)];
  const double Q22 = Z[(i + 1) * ySize + (j + 1)];

  // Interpolate in x for the "bottom" and "top" edges
  const double R1 = (1.0 - tx) * Q11 + tx * Q21;  // along j row
  const double R2 = (1.0 - tx) * Q12 + tx * Q22;  // along j+1 row

  // Interpolate those results along y
  return (1.0 - ty) * R1 + ty * R2;
}

}  // namespace utils

#endif