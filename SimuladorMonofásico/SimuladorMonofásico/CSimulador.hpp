#ifndef	CSIMULADOR_HPP
#define CSIMULADOR_HPP

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <fstream>

#include "CGas.hpp"
#include "CWell.hpp"
#include "CGrid.hpp"
#include "CProps.hpp"
#include "CMatrix.hpp"
#include "CFluido.hpp"
#include "CGnuplot.hpp"
#include "CLiquido.hpp"
#include "CReservoir.hpp"
#include "CDiscretization.hpp"

class CSimulador {
public:
	CSimulador();

	void run();
private:
	const double PI = 3.141592;
	CWell* well;
	CGrid* grid;
	CFluido* fluido;
	CReservoir* reservoir;
	CDiscretization* discretization;

	std::vector<std::vector<double>> Pressure;	/// todas as pressoes em todo o tempo
	std::vector<double> WellPressure;			/// todas as pressoes do poco

private:
	std::vector<double> calc_H(CGrid* grid, CProps* props_n, CProps* props_nu, double dt);
	std::vector<double> calc_Q(CGrid* grid, CWell* well, double time);
	std::vector<double> calc_X(double pw, std::vector<double> p);
	std::vector<std::vector<double>> calc_T(CGrid* grid, CProps* props_n);

	std::vector<std::vector<double>> calc_eta(CGrid* grid, CProps* props_nu, double dt);
	std::vector<std::vector<double>> calc_tau(CGrid* grid, CWell* well, CProps* props_nu, std::vector<double> p_nu, double pw_nu);

	void plot(std::vector<double> time, std::vector<double> Pw);

	bool isErrorNotAcceptable(double dt, std::vector<double> R, CGrid* grid, CProps* props_nu, double q, CDiscretization* discretization);

	void read_data_and_start_objects(std::string nameFile);
};
#endif