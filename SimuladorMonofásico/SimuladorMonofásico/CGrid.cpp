#include "CGrid.hpp"

CGrid::CGrid(CReservoir* reservoir, CDiscretization* discretization, CWell* well) {
	nr = discretization->get_nr();
	nz = discretization->get_nz();
	dtheta = reservoir->get_theta();
	createPosicoes(reservoir);
	createVolumes(reservoir, discretization->get_Ac());
	createPermeabilidades(reservoir, discretization);
	createFatorGeometricos(reservoir, discretization, well);
	createTime(discretization, well);
}

void CGrid::createPosicoes(CReservoir* reservoir) {
	alpha = pow(reservoir->get_re() / reservoir->get_rw(), 1.0 / nr);
	omega = log((alpha - 1) / log(alpha)) / log(alpha);

	/// -------------- raio --------------
	r.push_back(reservoir->get_rw() * log(alpha) / (1-(1/alpha)));
	rmh.push_back(reservoir->get_rw());
	for (int i = 1; i < nr; i++) {
		r.push_back(r[0] * pow(alpha, i));
		rph.push_back((r[i] - r[i-1]) / log(alpha));
		rmh.push_back(rph[i-1]);
	}
	rph.push_back(reservoir->get_re());

	/// ---------- profundidade ----------
	double dz = reservoir->get_h() / nz;
	for (int j = 0; j < nz; j++) {
		zmh.push_back(dz * j);
		zph.push_back(dz * ((float)j + 1.0));
		z.push_back((zmh[j] + zph[j]) / 2.0);
	}
}

void CGrid::createVolumes(CReservoir* reservoir, double Ac) {
	double vb;
	double dz = reservoir->get_h() / nz;
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nr; i++) {
			vb = 0.5 * (pow(rph[i], 2) - pow(rmh[i], 2)) * reservoir->get_theta() * dz;
			Vb.push_back( vb );
			Vb_ac.push_back(vb * Ac);
		}
	}
}

void CGrid::createPermeabilidades(CReservoir* reservoir, CDiscretization* discretization) {
	int nrs = discretization->get_nrs();

	double k0r = reservoir->get_k0r();
	double k0z = reservoir->get_k0z();

	double krs = k0r / (reservoir->get_S() / log(rph[nrs] / reservoir->get_rw()) + 1.0);
	double kzs = k0z / (reservoir->get_S() / log(rph[nrs] / reservoir->get_rw()) + 1.0);

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
			if (reservoir->get_infVol() && k == nz-1) {
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
	std::vector<double> partial = well->get_partial();
	for (int k = 0; k < nz; k++) {
		for (int i = 0; i < nr; i++) {
		/// fator geometrico em i - 1/2
		if (i == 0)
			Gimh.push_back(0.0);
		else
			Gimh.push_back(discretization->get_Bc() * rmh[i] * kimh[i + k * nr] * reservoir->get_theta() * (zph[k] - zmh[k]) / (r[i] - r[i - 1]));

		/// fator geometrico em i - 1/2
		if (i == nr-1)
			Giph.push_back(0.0);
		else
			Giph.push_back(discretization->get_Bc() * rph[i] * kiph[i + k * nr] * reservoir->get_theta() * (zph[k] - zmh[k]) / (r[i+1] - r[i]));

		/// fator geometrico em j - 1/2
		if (k == 0)
			Gjmh.push_back(0.0);
		else
			Gjmh.push_back(discretization->get_Bc() * (pow(rph[i], 2) - pow(rmh[i], 2)) * kjmh[i + k * nr] * reservoir->get_theta() / (2 * (z[k] - z[k - 1])));
			
		/// fator geometrico em j + 1/2
		if (k == nz-1)
			Gjph.push_back(0.0);
		else
			Gjph.push_back(discretization->get_Bc() * (pow(rph[i], 2) - pow(rmh[i], 2)) * kjph[i + k * nr] * reservoir->get_theta() / (2 * (z[k+1] - z[k])));
			
			
		/// fator geometrico do poco
		if (i == 0) 
			Gw.push_back(partial[k] * discretization->get_Bc() * rmh[0] * kimh[k * nr] * reservoir->get_theta() * (zph[k] - zmh[k]) / (r[0] - rmh[0]));
		}
	}
}

void CGrid::createTime(CDiscretization* discretization, CWell* well) {
	std::vector<double> tp = well->get_tp();
	double dtmin = discretization->get_dtmin();
	time.push_back(tp[0]);
	for (int i = 1; i < tp.size(); i++) {
		for (int t = 0; t < discretization->get_ntp(); t++) {
			time.push_back(tp[i-1] + pow(10, t * (log10(tp[i] - tp[i - 1]) - log10(dtmin)) / (discretization->get_ntp()-1.0))*dtmin);
		}
	}
}