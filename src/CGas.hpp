#ifndef CGAS_HPP
#define CGAS_HPP

#include<math.h>
#include <string>
#include<iostream>
#include "CFluido.hpp"

#include <iostream>

class CGas : public CFluido{
public:
	CGas(double _cf, double _p0, double _mu, double _T0, double _T, double _Tpc, double _Ppc, double _Ma):
		cf{ _cf }, p0{ _p0 }, mu{ _mu }, T0{ _T0 }, T{ _T }, Tpc{ _Tpc }, Ppc{ _Ppc }, Ma{ _Ma }{}
	CGas() {}

public:
	const double DELTA = 1.0e-5;
	double calc_b(double p);
	double calc_rho(double p);
	double calc_dbdp(double p);
	double calc_mu(double p);
	double calc_dmudp(double p);

	std::string get_type() { return type; }

private:

	double Z_KIAM(double p);
	double mu_LGE(double rho);

	std::string type = "gas";

	double cf{ 0.00215094 };	/// compressibilidade do fluido na condi??o inicial [cm^2/kgf]
	double p0{ 1.0335123 };		/// pressao padrao [kgf/cm^2]
	double mu{ 0.0262317 };		/// viscosidade na condicao inicial [cp]
	double T0{ 288.75 };		/// temperatura absoluta padrao [K]
	double T{ 353.15 };			/// temperatura do fluido no reservat?rio [K]
	double Tpc{ 216.32 };		/// temperatura pseudocritica [K]
	double Ppc{ 46.34 };		/// pressao pseudocritica [kgf/cm^2]
	double Ma{ 20.3 };			/// massa molecular aparente [kg/kg-mol]

	const double R = 0.08478;	///constante universal dos gases (ANP) [(kgf/cm^2)*(m^3)/(kg-mol*K)]
};
#endif