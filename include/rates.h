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
  const double m_lgEnergyMax = 23;
  const size_t m_lgEnergySize = 1000;

  std::vector<double> m_lgE;
  std::vector<double> m_beta;
};

class AbsorptionRate {
 public:
  AbsorptionRate(const std::string& filename);
  double get(double energy, double z = 0) const;

 private:
  const double m_zMin = 0;
  const double m_zMax = 10;
  const size_t m_zSize = 51;

  const double m_lgEnergyMin = 10;
  const double m_lgEnergyMax = 23;
  const size_t m_lgEnergySize = 1000;

  std::vector<double> m_z, m_lgE;
  std::vector<std::vector<double>> m_rate;
};

class ProductionRate {
 public:
  ProductionRate(const std::string& filename);
  double get(double E_sec, double E_pri, double z = 0) const;

 private:
  const double m_zMin = 0;
  const double m_zMax = 10;
  const size_t m_zSize = 51;

  const double m_lgxMin = -5;
  const double m_lgxMax = 0;
  const size_t m_xSize = 100;

  const double m_lgEnergyMin = 17;
  const double m_lgEnergyMax = 23;
  const size_t m_lgEnergySize = 1000;

  std::vector<double> m_z, m_lgx, m_lgE;
  std::vector<std::vector<double>> m_rate;
};

}  // namespace gnuprop

#endif  // GNUPROP_RATES_H