#ifndef CGRID_HPP
#define CGRID_HPP

#include  <vector>
#include <math.h>

#include "CReservoir.hpp"
#include "CDiscretization.hpp"
#include "CWell.hpp"

class CGrid {
public:
	CGrid(CReservoir* reservoir, CDiscretization* discretization, CWell* well);

private:
	int nr, nz;
	double alpha;
	double omega;
	double dtheta;

	std::vector<double> time;	/// vetor dos tempos

	std::vector<double> r;		/// posicao do centro do bloco em relacao a x
	std::vector<double> rph;	/// posicao em x + 1/2
	std::vector<double> rmh;	/// posicao em x - 1/2

	std::vector<double> z;		/// posicao do centro do bloco em relacao a z
	std::vector<double> zph;	/// posicao em z + 1/2
	std::vector<double> zmh;	/// posicao em z - 1/2

	std::vector<double> Vb;		/// volume 
	std::vector<double> Vb_ac;	/// volume corrigido

	std::vector<double> kr;		/// permeabilidade horizontal
	std::vector<double> kz;		/// permeabilidade vertical

	std::vector<double> kiph;	/// permeabilidade em i + 1/2
	std::vector<double> kimh;	/// permeabilidade em i - 1/2
	std::vector<double> kjph;	/// permeabilidade em j + 1/2
	std::vector<double> kjmh;	/// permeabilidade em j - 1/2

	std::vector<double> Giph;	/// fator geometrico da transmissibilidade em i + 1/2
	std::vector<double> Gimh;	/// fator geometrico da transmissibilidade em i - 1/2
	std::vector<double> Gjph;	/// fator geometrico da transmissibilidade em j + 1/2
	std::vector<double> Gjmh;	/// fator geometrico da transmissibilidade em j - 1/2

	std::vector<double> Gw;		/// fator geometrico para calculo da c.c. interna (poco)

private: /// private functions
	void createPosicoes(CReservoir* reservoir);
	void createVolumes(CReservoir* reservoir, double Ac);
	void createPermeabilidades(CReservoir* reservoir, CDiscretization* discretization);
	void createFatorGeometricos(CReservoir* reservoir, CDiscretization* discretization, CWell* well);
	void createTime(CDiscretization* discretization, CWell* well);

public:
	/// funcoes get
	double get_time( int i)	{ return time[i]; }

	double get_r( int i)	{ return r[i]; }
	double get_rph( int i)	{ return rph[i]; }
	double get_rmh( int i)	{ return rmh[i]; }

	double get_z( int i)	{ return z[i]; }
	double get_zph( int i)	{ return zph[i]; }
	double get_zmh( int i)	{ return zmh[i]; }

	double get_Vb( int i)	{ return Vb[i]; }
	double get_Vb_ac( int i){ return Vb_ac[i]; }

	double get_kr( int i)	{ return kr[i]; }
	double get_kz( int i)	{ return kz[i]; }

	double get_kiph( int i)	{ return kiph[i]; }
	double get_kimh( int i)	{ return kimh[i]; }
	double get_kjph( int i)	{ return kjph[i]; }
	double get_kjmh( int i)	{ return kjmh[i]; }

	double get_Giph( int i)	{ return Giph[i]; }
	double get_Gimh( int i)	{ return Gimh[i]; }
	double get_Gjph( int i)	{ return Gjph[i]; }
	double get_Gjmh( int i)	{ return Gjmh[i]; }

	double get_Gw( int i)	{ return Gw[i]; }

	/// get vetor completo
	std::vector<double> get_time() { return time; }

	std::vector<double> get_r() { return r; }
	std::vector<double> get_rph() { return rph; }
	std::vector<double> get_rmh() { return rmh; }

	std::vector<double> get_z() { return z; }
	std::vector<double> get_zph() { return zph; }
	std::vector<double> get_zmh() { return zmh; }

	std::vector<double> get_Vb() { return Vb; }
	std::vector<double> get_Vb_ac() { return Vb_ac; }

	std::vector<double> get_kr() { return kr; }
	std::vector<double> get_kz() { return kz; }

	std::vector<double> get_kiph() { return kiph; }
	std::vector<double> get_kimh() { return kimh; }
	std::vector<double> get_kjph() { return kjph; }
	std::vector<double> get_kjmh() { return kjmh; }

	std::vector<double> get_Giph() { return Giph; }
	std::vector<double> get_Gimh() { return Gimh; }
	std::vector<double> get_Gjph() { return Gjph; }
	std::vector<double> get_Gjmh() { return Gjmh; }

	std::vector<double> get_Gw() { return Gw; }


	int get_ntotal()				{ return nr * nz; }
	int get_nr()					{ return nr; }
	int get_nz()					{ return nz; }

	double get_alpha()				{ return alpha; }
	double get_omega()				{ return omega; }
	double get_dtheta()				{ return dtheta; }
};
#endif