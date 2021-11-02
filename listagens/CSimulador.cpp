#include "CSimulador.hpp"

CSimulador::CSimulador() {
	std::cout << "Nome do arquivo: ";
	std::string nameFile;
	nameFile = "input.dat";
	std::cout << nameFile << std::endl;
	//std::cin >> nameFile;
	read_data_and_start_objects(nameFile);
}

void CSimulador::run() {
	CProps* props_n  = new CProps(grid, fluido, reservoir, discretization);
	CProps* props_nu = new CProps(grid, fluido, reservoir, discretization);
	std::vector<double> time = grid->get_time();

	/// var de ajuda
	int nr = discretization->get_nr();
	int nz = discretization->get_nz();
	int iteracoes = 0;
	double erro_MB = 1;
	double erro_NR = 1;
	double dt;

	/// variaveis das pressoes - pressao no tempo anterior / pressao iteracao n / pressao iteracao n+1
	double pw_n		= reservoir->get_p_i();
	double pw_nu	= reservoir->get_p_i();
	std::vector<double> Pw(time.size());
	std::vector<double> p_n(nr * nz, reservoir->get_p_i());
	std::vector<double> p_nu(nr * nz, reservoir->get_p_i());

	/// comeco a preparar o loop
	std::vector<double> H;
	std::vector<double> Q;
	std::vector<double> X;
	std::vector<std::vector<double>> T;
	std::vector<double> R(nr*nz+1);

	std::vector<std::vector<double>> eta(nr*nz+1, std::vector<double>(nr*nz+1));
	std::vector<std::vector<double>> tau(nr * nz + 1, std::vector<double>(nr * nz + 1));
	std::vector<std::vector<double>> J(nr * nz + 1, std::vector<double>(nr * nz + 1));
	std::vector<double> dX;

	Pw[0] = reservoir->get_p_i();

	clock_t begin_time;
	/// loop do tempo
	for (unsigned int t = 1; t < time.size(); t++) {
		dt = time[t] - time[t - 1];

		props_n->update(pw_n, p_n);

		begin_time = clock();

		do {
			props_nu->update(pw_nu, p_nu);

			H = calc_H(grid, props_n, props_nu, dt);
			Q = calc_Q(grid, well, time[t]);
			X = calc_X(pw_nu, p_nu);
			T = calc_T(grid, props_nu);

			for (int i = 0; i < nr * nz + 1; i++) { /// a linha da multiplicacao
				R[i] = 0.0;
				for (int j = 0; j < nr * nz + 1; j++) /// os termos da multiplicacao
					R[i] += T[i][j] * X[j];
				R[i] += Q[i] - H[i];
			}


			/// matriz jacobiana
			eta = calc_eta(grid, props_nu, dt);
			tau = calc_tau(grid, well, props_nu, p_nu, pw_nu);

			for (int i = 0; i < nr * nz + 1; i++) /// a linha da multiplicacao
				for (int j = 0; j < nr * nz + 1; j++) /// os termos da multiplicacao
					J[i][j] = T[i][j] + tau[i][j] - eta[i][j];

			dX = CMatrix::LU_solver(J, R);
			//dX = CMatrix::GaussSolver(R, J); /// calculo a solucao do sistema linear

			pw_nu = pw_nu + dX[0];
			for (int i = 0; i < nr * nz; i++)
				p_nu[i] = p_nu[i] + dX[i+1];

			iteracoes++;
		} while (iteracoes < 10 && isErrorNotAcceptable(dt, R, grid, props_nu, Q[0], discretization));

		std::cout << "Tempo: " << time[t] << " - iteracoes: " << iteracoes << " - duracao: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
		p_n = p_nu;
		pw_n = pw_nu;
		iteracoes = 0;
		Pw[t] = pw_nu;
		Pressure.push_back(p_n);
	}
	CMatrix::mostrarMatriz(Pw, "Pw");
	plot(time, Pw);
}

std::vector<double> CSimulador::calc_H(CGrid* grid, CProps* props_n, CProps* props_nu, double dt) {
	std::vector<double> H(grid->get_ntotal()+1,0.0);
	for (int i = 0; i < grid->get_ntotal(); i++)
		H[i+1] = grid->get_Vb_ac(i) * (props_nu->get_b(i)* props_nu->get_phi(i) - props_n->get_b(i)* props_n->get_phi(i)) / dt;
	return H;
}

std::vector<double> CSimulador::calc_Q(CGrid* grid, CWell* well, double time) {
	std::vector<double> Q(grid->get_ntotal() + 1, 0.0);
	Q[0] = well->get_qsc(time) * (grid->get_dtheta() / (2 * PI));
	return Q;
}

std::vector<double> CSimulador::calc_X(double pw, std::vector<double> p) {
	std::vector<double> X(p.size() + 1);
	X[0] = pw;
	for (unsigned int i = 1; i < X.size(); i++)
		X[i] = p[i-1];
	return X;
}

std::vector<std::vector<double>> CSimulador::calc_T(CGrid* grid, CProps* props_n) {
	int nt = grid->get_ntotal()+1;
	int nr = grid->get_nr();
	int nz = grid->get_nz();
	int i;

	double Right;
	double Left;
	double Bottom;
	double Top;
	double Well;

	std::vector<std::vector<double>> T(nt, std::vector<double>(nt, 0.0));

	T[0].resize(nt, 0.0);
	for (int k = 1; k < nt; k++) { // rodo sem o poco
		i = k - 1;
		//T[k].resize(nt, 0.0);
		Right = 0.0;
		Left = 0.0;
		Bottom = 0.0;
		Top = 0.0;
		Well = 0.0;

		/// passo por todas as colunas da matriz
		if (i - nr >= 0) {
			Top = grid->get_Gjmh(i) * props_n->get_bjmh(i) / props_n->get_mujmh(i);
			T[k][k - nr] = Top;
		}
		if (i + nr < nt-1) {
			Bottom = grid->get_Gjph(i) * props_n->get_bjph(i) / props_n->get_mujph(i);
			T[k][k + nr] = Bottom;
		}
		if ((i - 1) % nr >= 0) {
			Left = grid->get_Gimh(i) * props_n->get_bimh(i) / props_n->get_muimh(i);
			T[k][k - 1] = Left;
		}
		if ((i + 1) % nr > 0) {
			Right = grid->get_Giph(i) * props_n->get_biph(i) / props_n->get_muiph(i);
			T[k][k + 1] = Right;
		}
		if (i % nr == 0) {
			Well = well->get_partial((int)i/nr) * grid->get_Gw((int)i/nr) * props_n->get_bw() / props_n->get_muw();
			T[k][0] = Well;
			T[0][k] = Well;
			T[0][0] -= Well;
		}
		T[k][k] = -Top - Bottom - Right - Left - Well;
	}
	return T;
}

std::vector<std::vector<double>> CSimulador::calc_eta(CGrid* grid, CProps* props_nu, double dt) {
	int nt = grid->get_ntotal() + 1;
	int nr = grid->get_nr();
	int nz = grid->get_nz();

	std::vector<std::vector<double>> eta(nt, std::vector<double> (nt, 0.0));
	eta[0].resize(nt, 0.0);
	for (int i = 1; i < nt; i++) {
		eta[i][i] = grid->get_Vb_ac(i-1)*(props_nu->get_b(i-1)*props_nu->get_dphidp(i-1) + props_nu->get_phi(i-1) * props_nu->get_dbdp(i-1)) / dt;
	}
	return eta;
}

std::vector<std::vector<double>> CSimulador::calc_tau(CGrid* grid, CWell* well, CProps* props_nu, std::vector<double> p_nu, double pw_nu) {
	int nt = grid->get_ntotal();
	int nr = grid->get_nr();
	int nz = grid->get_nz();

	double omega = grid->get_omega();

	std::vector<std::vector<double>> tau(nt + 1, std::vector<double>(nt+1, 0.0));

	/// bottom
	for (int i = 0; i < nt - nr; i++) {
		tau[i + 1][i + 1 + nr] = (p_nu[i] - p_nu[i + nr]) * grid->get_Gjph(i)
			* (props_nu->get_dbdp(i) * (props_nu->get_mu(i + nr) + props_nu->get_mu(i))
				- props_nu->get_dmudp(i) * (props_nu->get_b(i + nr) + props_nu->get_b(i)))
			/ pow(props_nu->get_mu(i) + props_nu->get_mu(i + nr), 2);
		tau[i + 1][i + 1] += tau[i + 1][i + 1 + nr];
	}

	/// top
	for (int i = nr; i < nt; i++) {
		tau[i + 1][i + 1 - nr] = (p_nu[i] - p_nu[i - nr]) * grid->get_Gjmh(i)
			* (props_nu->get_dbdp(i) * (props_nu->get_mu(i - nr) + props_nu->get_mu(i))
				- props_nu->get_dmudp(i) * (props_nu->get_b(i - nr) + props_nu->get_b(i)))
			/ pow(props_nu->get_mu(i) + props_nu->get_mu(i - nr), 2);
		tau[i + 1][i + 1] += tau[i + 1][i + 1 - nr];
	}

	/// left
	for (int i = 1; i < nt; i++) {
		tau[i+1][i] = grid->get_Gimh(i) * (1.0 - omega) * (p_nu[i - 1] - p_nu[i])
			* (props_nu->get_dbdp(i) / props_nu->get_muimh(i)
				- (props_nu->get_bimh(i) / pow(props_nu->get_muimh(i), 2)) * props_nu->get_dmudp(i - 1));
		tau[i+1][i+1] = tau[i+1][i];
	}

	/// right
	for (int i = 1; i < nt - 1; i++) {
		tau[i][i + 1] = grid->get_Giph(i) * omega * (p_nu[i + 1] - p_nu[i])
			* (props_nu->get_dbdp(i) / props_nu->get_muiph(i)
				- (props_nu->get_biph(i) / pow(props_nu->get_muiph(i), 2)) * props_nu->get_dmudp(i - 1));
		tau[i][i] = tau[i][i + 1];
	}

	/// well
	double temp_well;
	for (int i = 0; i < nz; i++) {
		tau[0][0] += (grid->get_Gw(i) / props_nu->get_muw()) * (props_nu->get_dbwdpw() - (props_nu->get_bw() / props_nu->get_muw()) * props_nu->get_dmuwdpw())*(p_nu[i*nr] - pw_nu) * well->get_partial(i);
		temp_well = (grid->get_Gw(i) / props_nu->get_muw()) * (props_nu->get_dbwdp1() - (props_nu->get_bw() / props_nu->get_muw()) * props_nu->get_dmuwdp1()) * (p_nu[i * nr] - pw_nu) * well->get_partial(i);
		tau[0][i * nr + 1] += temp_well;
		tau[0][0] -= temp_well;
		tau[i * nr + 1][i * nr + 1] -= temp_well;
		tau[i * nr + 1][0] += (grid->get_Gw(i) / props_nu->get_muw()) * (props_nu->get_dbwdp1() - (props_nu->get_bw() / props_nu->get_muw()) * props_nu->get_dmuwdp1()) * (p_nu[i * nr] - pw_nu) * well->get_partial(i);
	}
	return tau;
}

void CSimulador::plot(std::vector<double> time, std::vector<double> Pw) {
	std::string name = ("Pw_versus_time");

	std::ofstream outdata; //save data
	outdata.open((name + ".dat").c_str());
	outdata << "# time Temperature " << std::endl;
	for (int i = 0; i < Pw.size(); i++)
		outdata << time[i] << " " << Pw[i] << std::endl;

	CGnuplot::semilogx((name + ".dat").c_str(), "time", "Pw", (name + ".png").c_str());
}

bool CSimulador::isErrorNotAcceptable(double dt, std::vector<double> R, CGrid* grid, CProps* props_nu, double q, CDiscretization* discretization) {
	double sum_R_MB = R[0], sum_Divisor_MB = 0.0, NR = 0.0;
	for (int i = 0; i < R.size()-1; i++) {
		sum_R_MB += R[i+1];
		sum_Divisor_MB += (grid->get_Vb_ac(i) * props_nu->get_phi(i));
		NR += R[i + 1] / (grid->get_Vb_ac(i) * props_nu->get_phi(i));
	}
	double MB = dt * abs(sum_R_MB) / sum_Divisor_MB;

	double biggestNR = abs(R[0] / q) > NR ? abs(R[0] / q) : NR;
	return (MB > discretization->get_eps_MB() || biggestNR > discretization->get_eps_NR());
}

void CSimulador::read_data_and_start_objects(std::string nameFile) {
	std::ifstream file(nameFile);

	std::string text;
	std::getline(file, text);
	std::getline(file, text);
	std::getline(file, text);
	std::getline(file, text);
	std::getline(file, text);
	
	/// Pegando as variaveis do Poco
	std::getline(file, text); 

	std::vector<double> tp;
	std::vector<double> qsc;
	std::vector<double> dz;
	std::vector<double> partial;

	int periodos;
	file >> text; file >> text;
	periodos = std::stoi(text);

	/// pego os valores dos tps
	file >> text;
	for (int i = 0; i < periodos; i++) {
		file >> text;
		tp.push_back(std::stof(text));
	}
	std::getline(file, text);

	/// pego os valores dos qsc
	file >> text;
	for (int i = 0; i < periodos; i++) {
		file >> text;
		qsc.push_back(std::stof(text));
	}
	std::getline(file, text);

	/// pego os valores do h e partial
	file >> text; file >> text;
	periodos = std::stoi(text);
	std::getline(file, text);
	std::getline(file, text);

	for (int i = 0; i < periodos; i++) {
		file >> text;
		dz.push_back(std::stof(text));
		file >> text; file >> text;
		partial.push_back(std::stof(text));
	}

	file >> text; file >> text;
	double rw = std::stof(text);
	well = new CWell(tp, qsc, dz, partial, rw);

	///
	/// RESERVOIR
	///
	std::getline(file, text);	std::getline(file, text); std::getline(file, text);
	file >> text; file >> text;
	bool isLiquid = std::stoi(text);

	file >> text; file >> text;
	double re = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double theta = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double k0r = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double k0z = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double cphi = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double phi0 = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double p0 = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double p_i = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double S = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double Temperature = std::stof(text); std::getline(file, text);

	reservoir = new CReservoir(isLiquid, rw, re, dz, theta, k0r, k0z, cphi, phi0, p0, p_i, S, Temperature);

	///
	/// DISCRETIZACAO
	/// 
	std::getline(file, text); std::getline(file, text);
	int nz = dz.size();

	file >> text; file >> text;
	int nr = std::stoi(text); std::getline(file, text);

	file >> text; file >> text;
	int nrs = std::stoi(text); std::getline(file, text);

	file >> text; file >> text;
	int nt = std::stoi(text); std::getline(file, text);

	file >> text; file >> text;
	int ntp = std::stoi(text); std::getline(file, text);

	file >> text; file >> text;
	int max_iter = std::stoi(text); std::getline(file, text);

	file >> text; file >> text;
	double dtmin = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double eps_NR = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double eps_MB = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double Ac = std::stof(text); std::getline(file, text);

	file >> text; file >> text;
	double Bc = std::stof(text); std::getline(file, text);

	discretization = new CDiscretization(nz, nr, nrs, nt, ntp, max_iter, dtmin, eps_NR, eps_MB, Ac, Bc);

	/// 
	/// FLUID
	/// 
	if (reservoir->isLiquid())
		fluido = new CLiquido;
	else
		fluido = new CGas;

	grid = new CGrid(reservoir, discretization, well);
}
