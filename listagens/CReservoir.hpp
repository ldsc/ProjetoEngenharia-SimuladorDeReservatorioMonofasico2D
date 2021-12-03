#ifndef CRESERVOIR_HPP
#define CRESERVOIR_HPP

#include<string>
#include<vector>

class CReservoir {
private:
	bool _isLiquid{1};

	std::vector<double> dz{ 1, 1 };			/// altura do reservatorio
	double rw{ 0.09486 };		/// raio interno (poco)
	double re{ 3000.0 };			/// raio externo
	double theta{3.141592/6};	/// angulo estudado do reservatorio
	double k0r{500};				/// permeabilidade horizontal
	double k0z{ 100 };			/// permeabilidade vertical
	double cphi{1.0e-4};			/// compressibilidadae da formacao
	double phi0{ 0.2 };			/// porosidade inicial
	double p0{1.033512};			/// pressao de referencia
	double p_i{350.0};			/// pressoa inicial
	double s{ 0 };				/// fator de pelicula
	double temperature{353.15};	/// temperatura do reservatorio

	bool aquifer{false};		/// presenca de aquifero por fronteira de Neumann?
	bool infVol{false};			/// presenca de aquifero por volume infinito?

public:
	CReservoir() {}
	CReservoir( bool _isLiquid, double _rw, double _re, std::vector<double> _dz, double _theta, double _k0r, double _k0z, double _cphi, double _phi0, double _p0, double _p_i, double _S, double _Temperature) : _isLiquid{ _isLiquid }, 
		rw{ _rw }, re{ _re }, dz{ _dz }, theta{ _theta }, k0r{ _k0r }, k0z{ _k0z }, cphi{ _cphi }, phi0{ _phi0 }, p0{ _p0 }, p_i{ _p_i }, s{ _S }, temperature{ _Temperature } {}

	double calc_phi(double p);
	double calc_dphidp(double p);

	/// funcoes get
	bool isLiquid() { return _isLiquid; }
	double Rw() { return rw; }
	double Re() { return re; }
	std::vector<double> Dz() { return dz; }
	double Theta() { return theta; }
	double K0r() { return k0r; }
	double K0z() { return k0z; }
	double Cphi() { return cphi; }
	double Phi0() { return phi0; }
	double P0() { return p0; }
	double P_i() { return p_i; }
	double S(){ return s; }
	double Temperature() { return temperature; }
	bool Aquifer() { return aquifer; }
	bool InfVol() { return infVol; }
};
#endif