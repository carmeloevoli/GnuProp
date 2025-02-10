#ifndef GNUPROP_RATES_H
#define GNUPROP_RATES_H

#include <string>
#include <vector>

namespace gnuprop {

class ProtonLossRate {
 public:
  ProtonLossRate(const std::string& filename);

  double beta(double E, double z = 0) const;

 private:
  const double m_lgEnergyMin = 16;
  const double m_lgEnergyMax = 24;
  const size_t m_energySize = 1000;

  std::vector<double> m_lgE;

  std::vector<double> m_beta;
};

class PhotoPionProductionRate {
 public:
  PhotoPionProductionRate(const std::string& filename);

  // Compute beta
  double get(double E_sec, double E_p, double z = 0) const;

 private:
  const double m_zMin = 0;
  const double m_zMax = 10;
  const size_t m_zSize = 21;

  const double m_lgXMin = -5;
  const double m_lgXMax = 0;
  const size_t m_xSize = 60;

  const double m_lgEnergyMin = 17;
  const double m_lgEnergyMax = 23;
  const size_t m_energySize = 200;

  std::vector<double> m_z;
  std::vector<double> m_lgx;
  std::vector<double> m_lgE;

  std::vector<std::vector<double>> m_rate;
};

class GammaAbsorptionRate {
 public:
  GammaAbsorptionRate(const std::string& filename);

  // Compute beta
  double get(double E_gamma, double z = 0) const;

 private:
  const double m_zMin = 0;
  const double m_zMax = 10;
  const size_t m_zSize = 101;

  const double m_lgEnergyMin = 9;
  const double m_lgEnergyMax = 20;
  const size_t m_energySize = 100;

  std::vector<double> m_z;
  std::vector<double> m_lgE;

  std::vector<std::vector<double>> m_rate;
};

}  // namespace gnuprop

#endif  // GNUPROP_RATES_H