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

	std::vector<double> vb;		/// volume 
	std::vector<double> vb_ac;	/// volume corrigido

	std::vector<double> kr;		/// permeabilidade horizontal
	std::vector<double> kz;		/// permeabilidade vertical

	std::vector<double> kiph;	/// permeabilidade em i + 1/2
	std::vector<double> kimh;	/// permeabilidade em i - 1/2
	std::vector<double> kjph;	/// permeabilidade em j + 1/2
	std::vector<double> kjmh;	/// permeabilidade em j - 1/2

	std::vector<double> giph;	/// fator geometrico da transmissibilidade em i + 1/2
	std::vector<double> gimh;	/// fator geometrico da transmissibilidade em i - 1/2
	std::vector<double> gjph;	/// fator geometrico da transmissibilidade em j + 1/2
	std::vector<double> gjmh;	/// fator geometrico da transmissibilidade em j - 1/2

	std::vector<double> gw;		/// fator geometrico para calculo da c.c. interna (poco)

private: /// private functions
	void createPosicoes(CReservoir* reservoir);
	void createVolumes(CReservoir* reservoir, double Ac);
	void createPermeabilidades(CReservoir* reservoir, CDiscretization* discretization);
	void createFatorGeometricos(CReservoir* reservoir, CDiscretization* discretization, CWell* well);
	void createTime(CDiscretization* discretization, CWell* well);

public:
	/// funcoes get
	std::vector<double> Time() { return time; }
	double Time( int i)	{ return time[i]; }

	double R( int i)	{ return r[i]; }
	double Rph( int i)	{ return rph[i]; }
	double Rmh( int i)	{ return rmh[i]; }

	double Z( int i)	{ return z[i]; }
	double Zph( int i)	{ return zph[i]; }
	double Zmh( int i)	{ return zmh[i]; }

	double Vb( int i)	{ return vb[i]; }
	double Vb_ac( int i){ return vb_ac[i]; }

	double Kr( int i)	{ return kr[i]; }
	double Kz( int i)	{ return kz[i]; }

	double Kiph( int i)	{ return kiph[i]; }
	double Kimh( int i)	{ return kimh[i]; }
	double Kjph( int i)	{ return kjph[i]; }
	double Kjmh( int i)	{ return kjmh[i]; }

	double Giph( int i)	{ return giph[i]; }
	double Gimh( int i)	{ return gimh[i]; }
	double Gjph( int i)	{ return gjph[i]; }
	double Gjmh( int i)	{ return gjmh[i]; }
					 
	double Gw( int i)	{ return gw[i]; }

	int Ntotal()		{ return nr * nz; }
	int Nr()			{ return nr; }
	int Nz()			{ return nz; }
	int Nt()			{ return time.size(); }

	double Alpha()		{ return alpha; }
	double Omega()		{ return omega; }
	double DTheta()		{ return dtheta; }
};
#endif