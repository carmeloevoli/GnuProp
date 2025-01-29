#ifndef GNUPROP_INTERACTIONS_KELNERAHARONIAN2008_H
#define GNUPROP_INTERACTIONS_KELNERAHARONIAN2008_H

#include <string>
#include <vector>

#include "simprop/core/units.h"

namespace KelnerAharonian2008 {

class SecondarySpectrum {
 public:
  SecondarySpectrum(const std::string& filename);
  virtual ~SecondarySpectrum() = default;

 public:
  virtual double Phi(double eta, double x) const = 0;

 protected:
  double xMinus(double eta) const;
  double xPlus(double eta) const;
  double B(double rho) const;
  double s(double rho) const;
  double delta(double rho) const;

  virtual double psi(double rho) const = 0;
  virtual double xPrimeMinus(double eta) const = 0;
  virtual double xPrimePlus(double eta) const = 0;

 protected:
  static constexpr double m_r = SI::pionMassC2 / SI::protonMassC2;
  static constexpr double m_r2 = m_r * m_r;
  static constexpr double m_eta_0 = 2. * m_r + m_r2;  // Eq. 16

  std::vector<double> m_lnrho_table;
  std::vector<double> m_s_table;
  std::vector<double> m_delta_table;
  std::vector<double> m_lnB_table;
};

class AntiNuMuSpectrum final : public SecondarySpectrum {
 public:
  AntiNuMuSpectrum();
  double Phi(double eta, double x) const override;

 protected:
  double psi(double rho) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

class NuMuSpectrum final : public SecondarySpectrum {
 public:
  NuMuSpectrum();
  double Phi(double eta, double x) const override;

 protected:
  double psi(double rho) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

class AntiNuElectronSpectrum final : public SecondarySpectrum {
 public:
  AntiNuElectronSpectrum();
  double Phi(double eta, double x) const override;

 protected:
  double psi(double rho) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

class NuElectronSpectrum final : public SecondarySpectrum {
 public:
  NuElectronSpectrum();
  double Phi(double eta, double x) const override;

 protected:
  double psi(double rho) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

struct NeutrinoProductionSpectrum {
  NuMuSpectrum numu;
  AntiNuMuSpectrum antiNumu;
  NuElectronSpectrum nue;
  AntiNuElectronSpectrum antiNue;

  double nu_mu(double eta, double x) const { return numu.Phi(eta, x); }
  double nu_e(double eta, double x) const { return nue.Phi(eta, x); }
  double barnu_mu(double eta, double x) const { return antiNumu.Phi(eta, x); }
  double barnu_e(double eta, double x) const { return antiNue.Phi(eta, x); }

  double Phi(double eta, double x) const {
    return numu.Phi(eta, x) + antiNumu.Phi(eta, x) + nue.Phi(eta, x) + antiNue.Phi(eta, x);
  }
};

class GammaSpectrum final : public SecondarySpectrum {
 public:
  GammaSpectrum();
  double Phi(double eta, double x) const override;

 protected:
  double psi(double rho) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

class ElectronSpectrum final : public SecondarySpectrum {
 public:
  ElectronSpectrum();
  double Phi(double eta, double x) const override;

 protected:
  double psi(double rho) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

class PositronSpectrum final : public SecondarySpectrum {
 public:
  PositronSpectrum();
  double Phi(double eta, double x) const override;

 protected:
  double psi(double rho) const override;
  double xPrimeMinus(double eta) const override;
  double xPrimePlus(double eta) const override;
};

struct PairSpectrum {
  ElectronSpectrum eMinus;
  PositronSpectrum ePlus;

  double e_minus(double eta, double x) const { return eMinus.Phi(eta, x); }
  double e_plus(double eta, double x) const { return ePlus.Phi(eta, x); }

  double Phi(double eta, double x) const { return eMinus.Phi(eta, x) + ePlus.Phi(eta, x); }
};

}  // namespace KelnerAharonian2008

#endif  // GNUPROP_INTERACTIONS_KELNERAHARONIAN2008_H
