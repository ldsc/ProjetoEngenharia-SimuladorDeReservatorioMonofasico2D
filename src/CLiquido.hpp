#ifndef CLIQUIDO_HPP
#define CLIQUIDO_HPP

#include <string>
#include "CFluido.hpp"

class CLiquido : public CFluido{
public:
	CLiquido() {}
	CLiquido (float _cf, float _b0, float _p0, float _mu, float _cmu) :cf{ _cf }, b0{ _b0 }, p0{ _p0 }, mu{ _mu }, cmu{ _cmu } {}

	double calc_b(double p);
	double calc_dbdp(double p);

	double calc_mu(double p);
	double calc_dmudp(double p);

	std::string get_type() { return type; }

private:
	std::string type = "liquid";

	double cf{ 14.7e-5 };		/// compressibilidade do fluido [cm^2/kgf]
	double b0{ 1.0 };			/// inverso do fator volume formacao na pressao p0 [m^3 std / m^3]
	double p0{ 1.0335123 };		/// pressao de referencia [kgf/cm^2]
	double mu{ 1.0 };			/// viscosidade [cp]
	double cmu{ 0.0 };
};
#endif
