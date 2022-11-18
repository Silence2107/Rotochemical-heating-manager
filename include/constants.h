#ifndef CONSTANTS_H
#define CONSTANTS_H

/// @brief various constants useful in current research <para></para>
/// contains namespace scientific, namespace conversion, namespace apr4, namespace ist_ns
namespace constants
{
	/// @brief useful constants from math and physics<para></para>
	/// contains const double G,const double Pi
	namespace scientific
	{
		/// @brief Gravitational constant in GeV^-2 units
		const double G = 6.709E-39;
		/// @brief Nucleon mass in GeV units
		const double M_N = 0.93956563;
		/// @brief Value of pi
		const double Pi = 3.14159265359;
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
		/// @brief Conversion from GeV to MeV
		const double gev_over_mev = 1E3;
		/// @brief Conversion from GeV^-1 to s
		const double gev_s = 1.5192E24;
	}
	/// @brief APR4 notable constants
	namespace apr4
	{
		/// @brief baryonic density (in fm^-3) limits in APR4. _low and _upp represent limits of EoS itself
		/// while _core_limit and _crust_limit represent phase transition boundaries
		const double nbar_low = 6.023E-13, nbar_core_limit = 9E-2, nbar_crust_limit = 2.096E-2, nbar_upp = 1.89;
		/// @brief pressure (in dyne cm^-2) limits in APR4. _low and _upp represent limits of EoS itself
		const double pressure_low = 1.003E17, pressure_upp = 8.3308E36;
		/// @brief energy density (in g cm^-3) limits in APR4. _low and _upp represent limits of EoS itself <para></para>
		/// while _core_limit represents phase transition boundary
		const double edensity_low = 1.000E3, edensity_upp = 7.4456E15, edensity_core_limit = 1.5197E14;
	}
	/// @brief IST for NS notable constants
	namespace ist_ns
	{
		/// @brief baryonic density (in fm^-3) limits in APR4. _low and _upp represent limits of EoS itself
		/// while _core_limit and _crust_limit represent to-from crust boundaries
		const double nbar_low = 2.3739996827636742E-11, nbar_upp = 2.3189838273277710, nbar_crust_limit = 9.9798029952044190E-2, nbar_core_limit = 9.9999913289570197E-2;
		/// @brief pressure (in MeV fm^-3) limits in IST. _low and _upp represent limits of EoS
		const double pressure_low = 5.1065210580102853E-14, pressure_upp = 228788.58970172083;
		/// @brief energy density (in MeV fm^-3) limits in IST. _low and _upp represent limits of EoS itself
		/// while _core_limit represent core boundary
		const double edensity_low = 2.2134491254971723E-8, edensity_upp = 15892.136580408434, edensity_core_limit = 94.214131003471735;
	}
}

#endif