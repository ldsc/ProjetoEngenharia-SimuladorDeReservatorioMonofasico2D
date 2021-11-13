#include "CSimulador.hpp"

CSimulador::CSimulador() {
	std::cout << "Nome do arquivo: ";
	std::string nameFile;
	std::cin >> nameFile;

	//nameFile = "input.dat";
	std::cout << nameFile << std::endl;
	//std::cin >> nameFile;
	read_data_and_start_objects(nameFile);
}

void CSimulador::run() {
	CProps* props_n  = new CProps(grid, fluido, reservoir, discretization);
	CProps* props_nu = new CProps(grid, fluido, reservoir, discretization);

	/// var de ajuda
	int nr = discretization->Nr();
	int nz = discretization->Nz();
	int iteracoes = 0;
	double erro_MB = 1;
	double erro_NR = 1;
	double dt;

	/// variaveis das pressoes - pressao no tempo anterior / pressao iteracao n / pressao iteracao n+1
	double pw_n		= reservoir->P_i();
	double pw_nu	= reservoir->P_i();
	std::vector<double> Pw(grid->Nt());
	std::vector<double> p_n(nr * nz, reservoir->P_i());
	std::vector<double> p_nu(nr * nz, reservoir->P_i());

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

	Pw[0] = reservoir->P_i();

	clock_t begin_time;
	/// loop do tempo
	for (unsigned int t = 1; t < grid->Nt(); t++) {
		dt = grid->Time(t) - grid->Time(t - 1);

		props_n->Update(pw_n, p_n);

		begin_time = clock();

		do {
			props_nu->Update(pw_nu, p_nu);

			H = calc_H( props_n, props_nu, dt);
			Q = calc_Q( grid->Time(t));
			X = calc_X(pw_nu, p_nu);
			T = calc_T( props_nu);

			for (int i = 0; i < nr * nz + 1; i++) { /// a linha da multiplicacao
				R[i] = 0.0;
				for (int j = 0; j < nr * nz + 1; j++) /// os termos da multiplicacao
					R[i] += T[i][j] * X[j];
				R[i] += Q[i] - H[i];
			}


			/// matriz jacobiana
			eta = calc_eta( props_nu, dt);
			tau = calc_tau( props_nu, p_nu, pw_nu);

			for (int i = 0; i < nr * nz + 1; i++) /// a linha da multiplicacao
				for (int j = 0; j < nr * nz + 1; j++) /// os termos da multiplicacao
					J[i][j] = T[i][j] + tau[i][j] - eta[i][j];

			dX = CMatrix::LU_solver(J, R);
			//dX = CMatrix::GaussSolver(R, J); /// calculo a solucao do sistema linear

			pw_nu = pw_nu + dX[0];
			for (int i = 0; i < nr * nz; i++)
				p_nu[i] = p_nu[i] + dX[i+1];

			iteracoes++;
		} while (iteracoes < 10 && isErrorNotAcceptable(dt, R, props_nu, Q[0]));

		std::cout << "Tempo: " << grid->Time(t) << " - iteracoes: " << iteracoes << " - duracao: " << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
		p_n = p_nu;
		pw_n = pw_nu;
		iteracoes = 0;
		Pw[t] = pw_nu;
		Pressure.push_back(p_n);
	}
	CMatrix::mostrarMatriz(Pw, "Pw");
	plot(grid->Time(), Pw);
	plotSurface();
}

std::vector<double> CSimulador::calc_H(CProps* props_n, CProps* props_nu, double dt) {
	std::vector<double> H(grid->Ntotal()+1,0.0);
	for (int i = 0; i < grid->Ntotal(); i++)
		H[i+1] = grid->Vb_ac(i) * (props_nu->B(i)* props_nu->Phi(i) - props_n->B(i)* props_n->Phi(i)) / dt;
	return H;
}

std::vector<double> CSimulador::calc_Q( double time) {
	std::vector<double> Q(grid->Ntotal() + 1, 0.0);
	Q[0] = well->Qsc(time) * (grid->DTheta() / (2 * PI));
	return Q;
}

std::vector<double> CSimulador::calc_X(double pw, std::vector<double> p) {
	std::vector<double> X(p.size() + 1);
	X[0] = pw;
	for (unsigned int i = 1; i < X.size(); i++)
		X[i] = p[i-1];
	return X;
}

std::vector<std::vector<double>> CSimulador::calc_T(CProps* props_n) {
	int nt = grid->Ntotal()+1;
	int nr = grid->Nr();
	int nz = grid->Nz();
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
			Top = grid->Gjmh(i) * props_n->Bjmh(i) / props_n->Mujmh(i);
			T[k][k - nr] = Top;
		}
		if (i + nr < nt-1) {
			Bottom = grid->Gjph(i) * props_n->Bjph(i) / props_n->Mujph(i);
			T[k][k + nr] = Bottom;
		}
		if ((i - 1) % nr >= 0) {
			Left = grid->Gimh(i) * props_n->Bimh(i) / props_n->Muimh(i);
			T[k][k - 1] = Left;
		}
		if ((i + 1) % nr > 0) {
			Right = grid->Giph(i) * props_n->Biph(i) / props_n->Muiph(i);
			T[k][k + 1] = Right;
		}
		if (i % nr == 0) {
			Well = well->Partial((int)i/nr) * grid->Gw((int)i/nr) * props_n->Bw() / props_n->Muw();
			T[k][0] = Well;
			T[0][k] = Well;
			T[0][0] -= Well;
		}
		T[k][k] = -Top - Bottom - Right - Left - Well;
	}
	return T;
}

std::vector<std::vector<double>> CSimulador::calc_eta( CProps* props_nu, double dt) {
	int nt = grid->Ntotal() + 1;
	int nr = grid->Nr();
	int nz = grid->Nz();

	std::vector<std::vector<double>> eta(nt, std::vector<double> (nt, 0.0));
	eta[0].resize(nt, 0.0);
	for (int i = 1; i < nt; i++) {
		eta[i][i] = grid->Vb_ac(i-1)*(props_nu->B(i-1)*props_nu->Dphidp(i-1) + props_nu->Phi(i-1) * props_nu->Dbdp(i-1)) / dt;
	}
	return eta;
}

std::vector<std::vector<double>> CSimulador::calc_tau( CProps* props_nu, std::vector<double> p_nu, double pw_nu) {
	int nt = grid->Ntotal();
	int nr = grid->Nr();
	int nz = grid->Nz();

	double omega = grid->Omega();

	std::vector<std::vector<double>> tau(nt + 1, std::vector<double>(nt+1, 0.0));

	/// bottom
	for (int i = 0; i < nt - nr; i++) {
		tau[i + 1][i + 1 + nr] = (p_nu[i] - p_nu[i + nr]) * grid->Gjph(i)
			* (props_nu->Dbdp(i) * (props_nu->Mu(i + nr) + props_nu->Mu(i))
				- props_nu->Dmudp(i) * (props_nu->B(i + nr) + props_nu->B(i)))
			/ pow(props_nu->Mu(i) + props_nu->Mu(i + nr), 2);
		tau[i + 1][i + 1] += tau[i + 1][i + 1 + nr];
	}

	/// top
	for (int i = nr; i < nt; i++) {
		tau[i + 1][i + 1 - nr] = (p_nu[i] - p_nu[i - nr]) * grid->Gjmh(i)
			* (props_nu->Dbdp(i) * (props_nu->Mu(i - nr) + props_nu->Mu(i))
				- props_nu->Dmudp(i) * (props_nu->B(i - nr) + props_nu->B(i)))
			/ pow(props_nu->Mu(i) + props_nu->Mu(i - nr), 2);
		tau[i + 1][i + 1] += tau[i + 1][i + 1 - nr];
	}

	/// left
	for (int i = 1; i < nt; i++) {
		tau[i+1][i] = grid->Gimh(i) * (1.0 - omega) * (p_nu[i - 1] - p_nu[i])
			* (props_nu->Dbdp(i) / props_nu->Muimh(i)
				- (props_nu->Bimh(i) / pow(props_nu->Muimh(i), 2)) * props_nu->Dmudp(i - 1));
		tau[i+1][i+1] = tau[i+1][i];
	}

	/// right
	for (int i = 1; i < nt - 1; i++) {
		tau[i][i + 1] = grid->Giph(i) * omega * (p_nu[i + 1] - p_nu[i])
			* (props_nu->Dbdp(i) / props_nu->Muiph(i)
				- (props_nu->Biph(i) / pow(props_nu->Muiph(i), 2)) * props_nu->Dmudp(i - 1));
		tau[i][i] = tau[i][i + 1];
	}

	/// well
	double temp_well;
	for (int i = 0; i < nz; i++) {
		tau[0][0] += (grid->Gw(i) / props_nu->Muw()) * (props_nu->Dbwdpw() - (props_nu->Bw() / props_nu->Muw()) * props_nu->Dmuwdpw())*(p_nu[i*nr] - pw_nu) * well->Partial(i);
		temp_well = (grid->Gw(i) / props_nu->Muw()) * (props_nu->Dbwdp1() - (props_nu->Bw() / props_nu->Muw()) * props_nu->Dmuwdp1()) * (p_nu[i * nr] - pw_nu) * well->Partial(i);
		tau[0][i * nr + 1] += temp_well;
		tau[0][0] -= temp_well;
		tau[i * nr + 1][i * nr + 1] -= temp_well;
		tau[i * nr + 1][0] += (grid->Gw(i) / props_nu->Muw()) * (props_nu->Dbwdp1() - (props_nu->Bw() / props_nu->Muw()) * props_nu->Dmuwdp1()) * (p_nu[i * nr] - pw_nu) * well->Partial(i);
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

void CSimulador::plotSurface() {
	std::vector<double> pressure = Pressure[Pressure.size() - 1];

	std::string name = ("P_versus_xz");
	std::ofstream outdata; //save data
	outdata.open((name + ".dat").c_str());

	outdata << "# x z Pressure " << std::endl;
	for (int z = 0; z < grid->Nz(); z++)
		for (int r = 0; r < grid->Nr(); r++)
			outdata << grid->R(r) << " " << grid->Z(z) << " " << pressure[grid->Nr()*z+r] << std::endl;

	CGnuplot::surfacePlot((name + ".dat").c_str(), (name + ".png").c_str());
}

bool CSimulador::isErrorNotAcceptable(double dt, std::vector<double> R, CProps* props_nu, double q) {
	double sum_R_MB = R[0], sum_Divisor_MB = 0.0, NR = 0.0;
	for (int i = 0; i < R.size()-1; i++) {
		sum_R_MB += R[i+1];
		sum_Divisor_MB += (grid->Vb_ac(i) * props_nu->Phi(i));
		NR += R[i + 1] / (grid->Vb_ac(i) * props_nu->Phi(i));
	}
	double MB = dt * abs(sum_R_MB) / sum_Divisor_MB;

	double biggestNR = abs(R[0] / q) > NR ? abs(R[0] / q) : NR;
	return (MB > discretization->Eps_MB() || biggestNR > discretization->Eps_NR());
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
	if (reservoir->isLiquid()) {

		std::getline(file, text); std::getline(file, text);
		file >> text; file >> text;
		double cf = std::stof(text); std::getline(file, text);

		file >> text; file >> text;
		double b0 = std::stof(text); std::getline(file, text);

		file >> text; file >> text;
		double p0 = std::stof(text); std::getline(file, text);

		file >> text; file >> text;
		double mu = std::stof(text); std::getline(file, text);

		file >> text; file >> text;
		double cmu = std::stof(text); std::getline(file, text);
		fluido = new CLiquido(cf, b0, p0, mu, cmu);
	}
	else {
		std::getline(file, text); std::getline(file, text);
		file >> text; file >> text;
		double cf = std::stof(text); std::getline(file, text);

		file >> text; file >> text;
		double p0 = std::stof(text); std::getline(file, text);

		file >> text; file >> text;
		double mu = std::stof(text); std::getline(file, text);

		file >> text; file >> text;
		double T0 = std::stof(text); std::getline(file, text);

		file >> text; file >> text;
		double Tpc = std::stof(text); std::getline(file, text);

		file >> text; file >> text;
		double Ppc = std::stof(text); std::getline(file, text);

		file >> text; file >> text;
		double Ma = std::stof(text); std::getline(file, text);
		fluido = new CGas(cf, p0, mu, T0, Temperature, Tpc, Ppc, Ma);
	}

	grid = new CGrid(reservoir, discretization, well);
}
