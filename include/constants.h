#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "../include/auxiliaries.h"

/// @brief various constants useful in astrophysics
namespace constants
{
	/// @brief notable physics/math constants
	namespace scientific
	{
		/// @brief Gravitational constant in GeV^-2 units
		const double G = 6.709E-39;
		/// @brief Nucleon mass in GeV units
		const double M_N = 0.93956563;
		/// @brief Electron mass in GeV units
		const double M_e = 0.0005109989461;
		/// @brief Muon mass in GeV units
		const double M_mu = 0.1056583745;
		/// @brief Value of pi
		const double Pi = 3.14159265359;
		/// @brief Stefan-Boltzmann constant in natural units
		const double Sigma = Pi * Pi / 60.0;
		/// @brief Nuclear matter saturation density in GeV^3
		const double N_sat = 0.0012345;

	}
	/// @brief conversion table
	namespace conversion
	{
		/// @brief Conversion from g/cm^3 to GeV^4 (e.g. energy density) in natural units
		const double g_over_cm3_gev4 = 4.362E-18;
		/// @brief Conversion from dyne/cm^2 to GeV^4 (e.g. pressure) in natural units
		const double dyne_over_cm2_gev4 = 4.818E-39;
		/// @brief Conversion from Mev/fm3 to GeV^4 (e.g. pressure/energy density) in natural units
		const double mev_over_fm3_gev4 = 7.716E-6;
		/// @brief fm^3 GeV^3 in natural units
		const double fm3_gev3 = 1.296E2;
		/// @brief km GeV in natural units
		const double km_gev = 5.06E18;
		/// @brief Conversion from GeV to M_sol in natural units
		const double gev_over_msol = 8.951E-58;
		/// @brief Conversion from GeV to MeV in natural units
		const double gev_over_mev = 1E3;
		/// @brief Conversion from GeV^-1 to s in natural units
		const double gev_s = 1.5192E24;
		/// @brief Conversion from GeV to K in natural units
		const double gev_over_k = 1.1604E13;
		/// @brief Conversion from erg to GeV in natural units
		const double erg_over_gev = 6.2415E2;
		/// @brief Conversion from Myr to s in natural units
		const double myr_over_s = 3.1536E13;
		/// @brief Conversion from erg/(cm^3 s) to GeV^5 in natural units
		const double erg_over_cm3_s_gev5 = 3.17E-63;
	}
	/// @brief predefined species
	namespace species
	{
		/// @brief Neutron species
		const auxiliaries::phys::Species neutron(auxiliaries::phys::Species::ParticleType::kNeutron, auxiliaries::phys::Species::ParticleClassification::kBaryon);
		/// @brief Proton species
		const auxiliaries::phys::Species proton(auxiliaries::phys::Species::ParticleType::kProton, auxiliaries::phys::Species::ParticleClassification::kBaryon);
		/// @brief Electron species
		const auxiliaries::phys::Species electron(auxiliaries::phys::Species::ParticleType::kElectron, auxiliaries::phys::Species::ParticleClassification::kLepton);
		/// @brief Muon species
		const auxiliaries::phys::Species muon(auxiliaries::phys::Species::ParticleType::kMuon, auxiliaries::phys::Species::ParticleClassification::kLepton);
		/// @brief Tau species
		const auxiliaries::phys::Species tau(auxiliaries::phys::Species::ParticleType::kTau, auxiliaries::phys::Species::ParticleClassification::kLepton);
		/// @brief Pion species
		const auxiliaries::phys::Species pion(auxiliaries::phys::Species::ParticleType::kPion, auxiliaries::phys::Species::ParticleClassification::kMeson);
		/// @brief Kaon species
		const auxiliaries::phys::Species kaon(auxiliaries::phys::Species::ParticleType::kKaon, auxiliaries::phys::Species::ParticleClassification::kMeson);
		/// @brief Up quark species
		const auxiliaries::phys::Species uquark(auxiliaries::phys::Species::ParticleType::kUquark, auxiliaries::phys::Species::ParticleClassification::kQuark);
		/// @brief Down quark species
		const auxiliaries::phys::Species dquark(auxiliaries::phys::Species::ParticleType::kDquark, auxiliaries::phys::Species::ParticleClassification::kQuark);
		/// @brief Strange quark species
		const auxiliaries::phys::Species squark(auxiliaries::phys::Species::ParticleType::kSquark, auxiliaries::phys::Species::ParticleClassification::kQuark);
	}
}

#endif