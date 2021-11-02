#ifndef CRESERVOIR_HPP
#define CRESERVOIR_HPP

#include<string>
#include<vector>

class CReservoir {
private:
	bool _isLiquid{1};

	double rw{ 0.09486 };		/// raio interno (poco)
	double re{ 3000.0 };			/// raio externo
	std::vector<double> dz{ 1, 1 };			/// altura do reservatorio
	double theta{3.141592/6};	/// angulo estudado do reservatorio
	double k0r{500};				/// permeabilidade horizontal
	double k0z{ 100 };			/// permeabilidade vertical
	double cphi{1.0e-4};			/// compressibilidadae da formacao
	double phi0{ 0.2 };			/// porosidade inicial
	double p0{1.033512};			/// pressao de referencia
	double p_i{350.0};			/// pressoa inicial
	double S{ 0 };				/// fator de pelicula
	double Temperature{353.15};	/// temperatura do reservatorio

	bool aquifer{false};		/// presenca de aquifero por fronteira de Neumann?
	bool infVol{false};			/// presenca de aquifero por volume infinito?

public:
	CReservoir() {}
	CReservoir( bool _isLiquid, double _rw, double _re, std::vector<double> _dz, double _theta, double _k0r, double _k0z, double _cphi, double _phi0, double _p0, double _p_i, double _S, double _Temperature) : _isLiquid{ _isLiquid }, 
		rw{ _rw }, re{ _re }, dz{ _dz }, theta{ _theta }, k0r{ _k0r }, k0z{ _k0z }, cphi{ _cphi }, phi0{ _phi0 }, p0{ _p0 }, p_i{ _p_i }, S{ _S }, Temperature{ _Temperature } {}

	double calc_phi(double p);

	double calc_dphidp(double p);

	/// funcoes get
	bool isLiquid() { return _isLiquid; }
	double get_rw() { return rw; }
	double get_re() { return re; }
	//double get_h() { return h; }
	std::vector<double> get_dz() { return dz; }
	double get_theta() { return theta; }
	double get_k0r() { return k0r; }
	double get_k0z() { return k0z; }
	double get_cphi() { return cphi; }
	double get_phi0() { return phi0; }
	double get_p0() { return p0; }
	double get_p_i() { return p_i; }
	double get_S(){ return S; }
	double get_Temperature() { return Temperature; }
	bool get_aquifer() { return aquifer; }
	bool get_infVol() { return infVol; }

	void set_rw(double _rw) { rw = _rw; }
	void set_re(double _re) { re = _re; }
	//void set_h(double _h) { h = _h; }

};
#endif