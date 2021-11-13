#ifndef CDISCRETIZATION_HPP
#define CDISCRETIZATION_HPP

class CDiscretization {

private:
	int nz{ 2 };					/// qtd de volumes na altura
	int nr{ 50 };					/// qtd de volumes na largura
	int nrs{ 1 };					/// qtd de volumes na regiao danificada
	int nt{ 100 };					/// qtd de tempos
	int ntp{ 100 };					/// qtd de tempos
	int max_iter{ 24 };				/// numero maximo de iteracoes
	double dtmin{ 1.0 / 3600 };		/// passo de tempo minimo [h]
	double eps_NR{ 1.0e-6 };		/// tolerancia de convergencia dos residuos
	double eps_MB{ 1.0e-8 };		/// tolerancia de convergencia do balanco de materiais
	double Ac{ 24 };				/// constante de conversao de unidades acumulo (ANP)
	double Bc{ 0.0083621472 };		/// constante de conversao de unidades fluxo (ANP)

public:
	CDiscretization() {}
	CDiscretization(int _nz, int _nr, int _nrs, int _nt, int _ntp, int _max_iter, double _dtmin, double _eps_NR, double _eps_MB, double _Ac, double _Bc):
		nz{ _nz }, nr{ _nr }, nrs{ _nrs }, nt{ _nt }, ntp{ _ntp }, max_iter{ _max_iter }, dtmin{ _dtmin }, eps_NR{ _eps_NR }, eps_MB{ _eps_MB }, Ac{ _Ac }, Bc{ _Bc }{}
	/// funcoes get
	int get_nz()		{ return nz; }
	int get_nr()		{ return nr; }
	int get_nrs()		{ return nrs; }
	int get_nt()		{ return nt; }
	int get_ntp()		{ return ntp; }
	int get_max_iter()	{ return max_iter; }
	double get_dtmin()	{ return dtmin; }
	double get_eps_NR() { return eps_NR; }
	double get_eps_MB() { return eps_MB; }
	double get_Ac()		{ return Ac; }
	double get_Bc()		{ return Bc; }
};
#endif