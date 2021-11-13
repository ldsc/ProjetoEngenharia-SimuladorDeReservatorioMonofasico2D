#include "CGrid.hpp"

CGrid::CGrid(CReservoir* reservoir, CDiscretization* discretization, CWell* well) {
	nr = discretization->Nr();
	nz = discretization->Nz();
	dtheta = reservoir->Theta();
	createPosicoes(reservoir);
	createVolumes(reservoir, discretization->Ac());
	createPermeabilidades(reservoir, discretization);
	createFatorGeometricos(reservoir, discretization, well);
	createTime(discretization, well);
}

void CGrid::createPosicoes(CReservoir* reservoir) {
	alpha = pow(reservoir->Re() / reservoir->Rw(), 1.0 / nr);
	omega = log((alpha - 1) / log(alpha)) / log(alpha);

	/// -------------- raio --------------
	r.push_back(reservoir->Rw() * log(alpha) / (1-(1/alpha)));
	rmh.push_back(reservoir->Rw());
	for (int i = 1; i < nr; i++) {
		r.push_back(r[0] * pow(alpha, i));
		rph.push_back((r[i] - r[i-1]) / log(alpha));
		rmh.push_back(rph[i-1]);
	}
	rph.push_back(reservoir->Re());

	/// ---------- profundidade ----------
	std::vector<double> dz = reservoir->Dz();
	for (int j = 0; j < nz; j++) {
		zmh.push_back(0.0+(j==0?0:dz[j]+zmh[j-1]));
		zph.push_back(zmh[j]+dz[j]);
		z.push_back((zmh[j] + zph[j]) / 2.0);
	}
}

void CGrid::createVolumes(CReservoir* reservoir, double Ac) {
	double _vb;
	std::vector<double> dz = reservoir->Dz();
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nr; i++) {
			_vb = 0.5 * (pow(rph[i], 2) - pow(rmh[i], 2)) * reservoir->Theta() * dz[k];
			vb.push_back( _vb );
			vb_ac.push_back(_vb * Ac);
		}
	}
}

void CGrid::createPermeabilidades(CReservoir* reservoir, CDiscretization* discretization) {
	int nrs = discretization->Nrs();

	double k0r = reservoir->K0r();
	double k0z = reservoir->K0z();

	double krs = k0r / (reservoir->S() / log(rph[nrs] / reservoir->Rw()) + 1.0);
	double kzs = k0z / (reservoir->S() / log(rph[nrs] / reservoir->Rw()) + 1.0);

	/// nos centros dos volumes
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nr; i++) {
			/// regiao danificada
			if (i <= nrs) {
				kr.push_back(k0r + krs);
				kz.push_back(k0z + kzs);
			}
			/// regiao normal
			else {
				kr.push_back(k0r);
				kz.push_back(k0z);
			}
			/// volume infinito
			if (reservoir->InfVol() && k == nz-1) {
				kr[i + k * nr] = kr[i + k * nr] * 1.0e4; /// multiplicador arbitrario, valor deve ser grande
				kz[i + k * nr] = kz[i + k * nr] * 1.0e4; /// multiplicador arbitrario, valor deve ser grande
			}

			/// permeabilidades homogeneas entre as camadas
			kiph.push_back(k0r);
			kimh.push_back(k0r);

			kjph.push_back(k0z);
			kjmh.push_back(k0z);
		}
	}
}

void CGrid::createFatorGeometricos(CReservoir* reservoir, CDiscretization* discretization, CWell* well) {
	double fatorGeometrico;
	std::vector<double> partial = well->Partial();
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nr; i++) {
		/// fator geometrico em i - 1/2
		if (i == 0)
			gimh.push_back(0.0);
		else
			gimh.push_back(discretization->Bc() * rmh[i] * kimh[i + k * nr] * reservoir->Theta() * (zph[k] - zmh[k]) / (r[i] - r[i - 1]));

		/// fator geometrico em i - 1/2
		if (i == nr-1)
			giph.push_back(0.0);
		else
			giph.push_back(discretization->Bc() * rph[i] * kiph[i + k * nr] * reservoir->Theta() * (zph[k] - zmh[k]) / (r[i+1] - r[i]));

		/// fator geometrico em j - 1/2
		if (k == 0)
			gjmh.push_back(0.0);
		else
			gjmh.push_back(discretization->Bc() * (pow(rph[i], 2) - pow(rmh[i], 2)) * kjmh[i + k * nr] * reservoir->Theta() / (2 * (z[k] - z[k - 1])));
			
		/// fator geometrico em j + 1/2
		if (k == nz-1)
			gjph.push_back(0.0);
		else
			gjph.push_back(discretization->Bc() * (pow(rph[i], 2) - pow(rmh[i], 2)) * kjph[i + k * nr] * reservoir->Theta() / (2 * (z[k+1] - z[k])));
			
			
		/// fator geometrico do poco
		if (i == 0) 
			gw.push_back(partial[k] * discretization->Bc() * rmh[0] * kimh[k * nr] * reservoir->Theta() * (zph[k] - zmh[k]) / (r[0] - rmh[0]));
		}
	}
}

void CGrid::createTime(CDiscretization* discretization, CWell* well) {
	std::vector<double> tp = well->Tp();
	double dtmin = discretization->DTmin();
	time.push_back(tp[0]);
	for (int i = 1; i < tp.size(); i++) {
		for (int t = 0; t < discretization->Ntp(); t++) {
			time.push_back(tp[i-1] + pow(10, t * (log10(tp[i] - tp[i - 1]) - log10(dtmin)) / (discretization->Ntp()-1.0))*dtmin);
		}
	}
}