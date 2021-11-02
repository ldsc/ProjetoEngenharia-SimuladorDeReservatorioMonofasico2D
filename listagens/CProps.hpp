#ifndef CPROPS_HPP
#define CPROPS_HPP

#include <vector>
#include "CGrid.hpp"
#include "CFluido.hpp"
#include "CReservoir.hpp"
#include "CDiscretization.hpp"

class CProps {
public:
	CProps(CGrid* _grid, CFluido* _fluido, CReservoir* _reservoir, CDiscretization* _discretization);

private:
	CGrid* grid;
	CFluido* fluido;
	CReservoir* reservoir;
	CDiscretization* discretization;

	int nr;
	int nz;

	std::vector<double> b;
	std::vector<double> dbdp;

	std::vector<double> biph;
	std::vector<double> bimh;
	std::vector<double> bjph;
	std::vector<double> bjmh;


	std::vector<double> mu;
	std::vector<double> dmudp;

	std::vector<double> muiph;
	std::vector<double> muimh;
	std::vector<double> mujph;
	std::vector<double> mujmh;


	std::vector<double> phi;
	std::vector<double> dphidp;

	double bw = .0;
	double dbwdp1 = .0;
	double dbwdpw = .0;
	double muw = .0;
	double dmuwdp1 = .0;
	double dmuwdpw = .0;

	void update_b(std::vector<double> pressure);
	void update_mu(std::vector<double> pressure);
	void update_well(double pw);

public:
	void update(double pw, std::vector<double> pressure);

	/// funcoes get do vetor completo
	std::vector<double> get_phi()		{ return phi; }
	std::vector<double> get_dphidp()	{ return dphidp; }
	std::vector<double> get_b()			{ return b; }
	std::vector<double> get_dbdp()		{ return dbdp; }

	std::vector<double> get_biph()		{ return biph; }
	std::vector<double> get_bimh()		{ return bimh; }
	std::vector<double> get_bjph()		{ return bjph; }
	std::vector<double> get_bjmh()		{ return bjmh; }


	std::vector<double> get_mu()		{ return mu; }
	std::vector<double> get_dmudp()		{ return dmudp; }

	std::vector<double> get_muiph()		{ return mujph; }
	std::vector<double> get_muimh()		{ return muimh; }
	std::vector<double> get_mujph()		{ return mujph; }
	std::vector<double> get_mujmh()		{ return mujmh; }

	/// get com as posicoes desejadas
	double get_phi(int i)		{ return phi[i]; }
	double get_dphidp(int i)	{ return dphidp[i]; }


	double get_b(int i)			{ return b[i]; }
	double get_dbdp(int i)		{ return dbdp[i]; }

	double get_biph(int i)		{ return biph[i]; }
	double get_bimh(int i)		{ return bimh[i]; }
	double get_bjph(int i)		{ return bjph[i]; }
	double get_bjmh(int i)		{ return bjmh[i]; }


	double get_mu(int i)		{ return mu[i]; }
	double get_dmudp(int i)		{ return dmudp[i]; }

	double get_muiph(int i)		{ return mujph[i]; }
	double get_muimh(int i)		{ return muimh[i]; }
	double get_mujph(int i)		{ return mujph[i]; }
	double get_mujmh(int i)		{ return mujmh[i]; }

	/// props do poco
	double get_bw()					{ return bw; }
	double get_dbwdp1()				{ return dbwdp1; }
	double get_dbwdpw()				{ return dbwdpw; }
	double get_muw()				{ return muw; }
	double get_dmuwdp1()			{ return dmuwdp1; }
	double get_dmuwdpw()			{ return dmuwdpw; }
};
#endif