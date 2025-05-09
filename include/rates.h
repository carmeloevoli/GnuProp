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
  double m_lgEnergyMin = 0;
  double m_lgEnergyMax = 0;
  size_t m_lgEnergySize = 0;

  std::vector<double> m_lgE;
  std::vector<double> m_beta;
};

class AbsorptionRate {
 public:
  AbsorptionRate(const std::string& filename);
  double get(double energy, double z = 0) const;

 private:
  double m_zMin = 0;
  double m_zMax = 0;
  size_t m_zSize = 0;

  double m_lgEnergyMin = 0;
  double m_lgEnergyMax = 0;
  size_t m_lgEnergySize = 0;

  std::vector<double> m_z, m_lgE;
  std::vector<std::vector<double>> m_rate;
};
;

class ProductionRate {
 public:
  ProductionRate(const std::string& filename);
  void setLimits(std::vector<double>& data);
  double get(double x, double E_pri, double z = 0) const;

 private:
  double m_zMin = 0;
  double m_zMax = 0;
  size_t m_zSize = 0;

  double m_lgXMin = 0;
  double m_lgXMax = 0;
  size_t m_lgXSize = 0;

  double m_lgEnergyMin = 0;
  double m_lgEnergyMax = 0;
  size_t m_lgEnergySize = 0;

  std::vector<double> m_z, m_lgX, m_lgE;
  std::vector<std::vector<double>> m_rate;
};

}  // namespace gnuprop

#endif  // GNUPROP_RATES_H