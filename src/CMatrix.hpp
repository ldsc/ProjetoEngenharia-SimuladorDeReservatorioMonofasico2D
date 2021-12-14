#ifndef CMATRIX_HPP
#define CMATRIX_HPP

#include<vector>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>

class CMatrix {
	/// esta classe nao oferece pilulas azuis nem vermelhas.
public:
	static double multp_vector(std::vector<double> A, std::vector<double> B);

	static std::vector<std::vector<double>> multp_matrix(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B);

	static std::vector<double> Moore_Penrose(std::vector<double> R, std::vector<std::vector<double>> J);
	static std::vector<double> GaussSolver(std::vector<double> R, std::vector<std::vector<double>> J);

	static std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> matrix);
	static double determinanteSparse(std::vector<std::vector<double>> matriz);
	static std::vector<double> LU_solver(std::vector<std::vector<double>> J, std::vector<double> R);
	static double maxAbs(std::vector<double> R);
	static std::vector<std::vector<double>> rand(int n);
	void static mostrarMatriz(std::vector<double> vector, std::string txt, bool stop = 0);
	static void mostrarMatriz(std::vector<std::vector<double>> matriz, std::string txt, bool stop = 0);
}; 
#endif