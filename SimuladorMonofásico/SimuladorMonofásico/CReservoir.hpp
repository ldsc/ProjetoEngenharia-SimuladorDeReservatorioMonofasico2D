#ifndef CRESERVOIR_HPP
#define CRESERVOIR_HPP

#include<string>

class CReservoir {
private:
	bool _isLiquid = 1;

	double rw{ 0.09486 };		/// raio interno (poco)
	double re{ 3000.0 };			/// raio externo
	double h{ 2.0 };			/// altura do reservatorio
	double theta{3.141592/6};	/// angulo estudado do reservatorio
	double k0r{500};				/// permeabilidade horizontal
	double k0z{ 100 };			/// permeabilidade vertical
	double cphi{1.0e-4};			/// compressibilidadae da formacao
	double phi0{ 0.2 };			/// porosidade inicial
	double p0{1.033512};			/// pressao de referencia
	double p_i{350.0};			/// pressoa inicial
	double S{ 0 };				/// fator de pelicula
	double Temperature{353.15};	/// temperatura do reservatorio
	double beta{ 3.1831588e+9 };	/// beta de Forchheimer

	bool aquifer{false};		/// presenca de aquifero por fronteira de Neumann?
	bool infVol{false};			/// presenca de aquifero por volume infinito?

public:
	double calc_phi(double p);

	double calc_dphidp(double p);

	/// funcoes get
	bool isLiquid() { return _isLiquid; }
	double get_rw() { return rw; }
	double get_re() { return re; }
	double get_h() { return h; }
	double get_theta() { return theta; }
	double get_k0r() { return k0r; }
	double get_k0z() { return k0z; }
	double get_cphi() { return cphi; }
	double get_phi0() { return phi0; }
	double get_p0() { return p0; }
	double get_p_i() { return p_i; }
	double get_S(){ return S; }
	double get_Temperature() { return Temperature; }
	double get_beta() { return beta; }
	bool get_aquifer() { return aquifer; }
	bool get_infVol() { return infVol; }

	void set_rw(double _rw) { rw = _rw; }
	void set_re(double _re) { re = _re; }
	void set_h(double _h) { h = _h; }

};
#endif