#ifndef INCLUDE_SIMPROP_UNITS_H
#define INCLUDE_SIMPROP_UNITS_H

#include <cmath>

#define pow2(A) ((A) * (A))
#define pow3(A) ((A) * (A) * (A))
#define pow4(A) ((A) * (A) * (A) * (A))

namespace SI {

// SI  base units
static constexpr double meter = 1;
static constexpr double second = 1;
static constexpr double metre = 1;
static constexpr double kilogram = 1;
static constexpr double ampere = 1;
static constexpr double kelvin = 1;
static constexpr double teradian = 1;

// SI derived units
static constexpr double radian = meter / meter;
static constexpr double hertz = 1 / second;
static constexpr double joule = kilogram * pow2(meter / second);
static constexpr double watt = joule / second;
static constexpr double steradian = pow2(meter / meter);

// prefixes
static constexpr double exa = 1e18;
static constexpr double peta = 1e15;
static constexpr double tera = 1e12;
static constexpr double giga = 1e9;
static constexpr double mega = 1e6;
static constexpr double kilo = 1e3;
static constexpr double milli = 1e-3;
static constexpr double micro = 1e-6;
static constexpr double nano = 1e-9;
static constexpr double pico = 1e-12;
static constexpr double femto = 1e-15;

// TIME UNITS
static constexpr double year = 3.154e7 * second;
static constexpr double kiloyear = kilo * year;
static constexpr double megayear = mega * year;
static constexpr double gigayear = giga * year;

// LENGTH UNITS
static constexpr double micron = micro * meter;
static constexpr double centimeter = 1e-2 * meter;
static constexpr double parsec = 3.086e16 * meter;
static constexpr double kiloparsec = kilo * parsec;
static constexpr double megaparsec = mega * parsec;
static constexpr double gigaparsec = giga * parsec;

// MASS UNITS
static constexpr double gram = 1e-3 * kilogram;

// ENERGY UNITS
static constexpr double electronvolt = 1.60217657e-19 * joule;
static constexpr double kiloelectronvolt = kilo * electronvolt;
static constexpr double megaelectronvolt = mega * electronvolt;
static constexpr double gigaelectronvolt = giga * electronvolt;
static constexpr double teraelectronvolt = tera * electronvolt;
static constexpr double petaelectronvolt = peta * electronvolt;
static constexpr double exaelectronvolt = exa * electronvolt;
static constexpr double erg = 1e-7 * joule;

// ABBREVIATION
static constexpr double sec = second;
static constexpr double kyr = kiloyear;
static constexpr double Myr = megayear;
static constexpr double Gyr = gigayear;
static constexpr double kpc = kiloparsec;
static constexpr double Mpc = megaparsec;
static constexpr double Mpc3 = Mpc * Mpc * Mpc;
static constexpr double eV = electronvolt;
static constexpr double keV = kiloelectronvolt;
static constexpr double MeV = megaelectronvolt;
static constexpr double GeV = gigaelectronvolt;
static constexpr double GeV2 = GeV * GeV;
static constexpr double TeV = teraelectronvolt;
static constexpr double PeV = petaelectronvolt;
static constexpr double cm2 = centimeter * centimeter;
static constexpr double cm3 = centimeter * centimeter * centimeter;
static constexpr double m2 = meter * meter;
static constexpr double m3 = meter * meter * meter;
static constexpr double km = kilo * meter;
static constexpr double Hz = hertz;
static constexpr double GHz = giga * hertz;
static constexpr double nW = nano * watt;
static constexpr double sr = steradian;
static constexpr double K = kelvin;

// PHYSICAL CONSTANTS
static constexpr double cLight = 2.99792458e8 * meter / second;
static constexpr double cLight2 = pow2(cLight);
static constexpr double cOver4pi = cLight / 4. / M_PI;
static constexpr double protonMass = 1.67262158e-24 * gram;
static constexpr double protonMassC = protonMass * cLight;
static constexpr double protonMassC2 = protonMass * cLight2;
static constexpr double neutronMass = 1.67492735e-24 * gram;
static constexpr double neutronMassC2 = neutronMass * cLight2;
static constexpr double electronMass = 9.10938291e-28 * gram;
static constexpr double electronMassC2 = electronMass * cLight2;
static constexpr double pionMassC2 = 134.9768 * MeV;
static constexpr double muonMassC2 = 105.6583755 * MeV;
static constexpr double sunMass = 1.989e33 * gram;
static constexpr double hPlanck = 6.62607015e-34 * joule * second;
static constexpr double hbarPlanck = hPlanck / 2. / M_PI;
static constexpr double hbarC = hbarPlanck * cLight;
static constexpr double kBoltzmann = 1.3806488e-23 * joule / kelvin;
static constexpr double electronRadius = 2.8179403227e-15 * meter;
static constexpr double barn = 1e-28 * m2;
static constexpr double mbarn = 1e-3 * barn;
static constexpr double alpha = 7.297352e-3;
static constexpr double sigmaTh = 6.6524e-25 * pow2(centimeter);

// COSMOLOGY
static constexpr double hLittle = 0.7;
static constexpr double H0 = 70. * hLittle * SI::km / SI::sec / SI::Mpc;
static constexpr double OmegaB = 0.;
static constexpr double OmegaC = 0.3;
static constexpr double OmegaM = OmegaC;
static constexpr double OmegaL = 1. - OmegaM;
static constexpr double Dh = SI::cLight / H0;

} // namespace SI

#endif // INCLUDE_SIMPROP_UNITS_H