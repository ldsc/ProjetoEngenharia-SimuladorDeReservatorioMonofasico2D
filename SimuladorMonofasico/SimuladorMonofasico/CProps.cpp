#include "CProps.hpp"

CProps::CProps(CGrid* _grid, CFluido* _fluido, CReservoir* _reservoir, CDiscretization* _discretization) {
	grid = _grid ; 
	fluido = _fluido;
	reservoir = _reservoir;
	discretization = _discretization;

	nr = discretization->Nr();
	nz = discretization->Nz();

	b.resize(nr*nz);
	dbdp.resize(nr * nz);

	biph.resize(nr * nz);
	bimh.resize(nr * nz);
	bjph.resize(nr * nz);
	bjmh.resize(nr * nz);

	mu.resize(nr * nz);
	dmudp.resize(nr * nz);

	muiph.resize(nr * nz);
	muimh.resize(nr * nz);
	mujph.resize(nr * nz);
	mujmh.resize(nr * nz);

	phi.resize(nr * nz);
	dphidp.resize(nr * nz);
}

void CProps::Update(double pw, std::vector<double> pressure) {
	for (int i = 0; i < nr*nz; i++) {
		b[i] = fluido->calc_b(pressure[i]);
		dbdp[i] = fluido->calc_dbdp(pressure[i]);

		mu[i] = fluido->calc_mu(pressure[i]);
		dmudp[i] = fluido->calc_dmudp(pressure[i]);

		phi[i] = reservoir->calc_phi(pressure[i]);
		dphidp[i] = reservoir->calc_dphidp(pressure[i]);
	}

	update_b(pressure);
	update_mu(pressure);
	update_well(pw);
}

void CProps::update_well(double pw) {
	bw = fluido->calc_b(pw);
	dbwdp1 = 0;
	dbwdpw = fluido->calc_dbdp(pw);
	muw = fluido->calc_mu(pw);
	dmuwdp1 = 0;
	dmuwdpw = fluido->calc_dmudp(pw);
}

void CProps::update_b(std::vector<double> pressure) {
	double omega = grid->Omega();
	int camada = -1;

	for (int i = 0; i < nr*nz; i++) {
		/// guardo a camada para uso futuro
		if (i % nr == 0) camada++;

		/// fronteira externa radial
		if ((i + 1) % nr == 0)
			biph[i] = b[i];
		else
			biph[i] = (1.0 - omega) * b[i] + omega * b[i + 1];

		/// fronteira interna radial
		if (i % nr == 0)
			bimh[i] = b[i];
		else
			bimh[i] = biph[i - 1];

		/// fronteira inferior vertical
		if (i >= (nz-1)*nr)
			bjph[i] = b[i];
		else
			bjph[i] = (b[i] * grid->Z(camada) + b[i+nr] * grid->Z(camada+1)) / (grid->Z(camada) + grid->Z(camada+1));
		/// fronteira superior vertical
		if (i < nr)
			bjmh[i] = b[i];
		else
			bjmh[i] = bjph[i-nr];
	}
}

void CProps::update_mu(std::vector<double> pressure) {
	double omega = grid->Omega();
	int camada = -1;

	for (int i = 0; i < nr * nz; i++) {
		/// guardo a camada para uso futuro
		if (i % nr == 0) camada++;

		/// fronteira externa radial
		if ((i + 1) % nr == 0)
			muiph[i] = mu[i];
		else
			muiph[i] = (1.0 - omega) * mu[i] + omega * mu[i + 1];

		/// fronteira interna radial
		if (i % nr == 0)
			muimh[i] = mu[i];
		else
			muimh[i] = muiph[i - 1];

		/// fronteira inferior vertical
		if (i >= (nz - 1) * nr)
			mujph[i] = mu[i];
		else
			mujph[i] = (mu[i] * grid->Z(camada) + mu[i + nr] * grid->Z(camada + 1)) / (grid->Z(camada) + grid->Z(camada + 1));
		/// fronteira superior vertical
		if (i < nr)
			mujmh[i] = mu[i];
		else
			mujmh[i] = mujph[i - nr];
	}
}