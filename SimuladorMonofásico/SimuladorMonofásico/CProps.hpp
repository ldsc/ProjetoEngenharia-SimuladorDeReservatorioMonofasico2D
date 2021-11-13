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
	void Update(double pw, std::vector<double> pressure);

	/// get com as posicoes desejadas
	double Phi(int i)		{ return phi[i]; }
	double Dphidp(int i)	{ return dphidp[i]; }


	double B(int i)			{ return b[i]; }
	double Dbdp(int i)		{ return dbdp[i]; }
	double Biph(int i)		{ return biph[i]; }
	double Bimh(int i)		{ return bimh[i]; }
	double Bjph(int i)		{ return bjph[i]; }
	double Bjmh(int i)		{ return bjmh[i]; }

	double Mu(int i)		{ return mu[i]; }
	double Dmudp(int i)		{ return dmudp[i]; }
	double Muiph(int i)		{ return mujph[i]; }
	double Muimh(int i)		{ return muimh[i]; }
	double Mujph(int i)		{ return mujph[i]; }
	double Mujmh(int i)		{ return mujmh[i]; }

	/// props do poco
	double Bw()					{ return bw; }
	double Dbwdp1()				{ return dbwdp1; }
	double Dbwdpw()				{ return dbwdpw; }
	double Muw()				{ return muw; }
	double Dmuwdp1()			{ return dmuwdp1; }
	double Dmuwdpw()			{ return dmuwdpw; }
};
#endif