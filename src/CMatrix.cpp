#include "CMatrix.hpp"

double CMatrix::multp_vector(std::vector<double> A, std::vector<double> B) {
	double result=0.0;
	for (int i = 0; i < A.size(); i++)
		result += A[i] * B[i];
	return result;
}

std::vector<std::vector<double>> CMatrix::multp_matrix(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B) {
	std::vector<std::vector<double>> result;
	result.resize(A.size());
	for (int i = 0; i < A.size(); i++)
		result[i].resize(B.size());

	double temp_value;
	for (int i = 0; i < A.size(); i++) {
		for (int k = 0; k < B[0].size(); k++) {
			temp_value = 0;
			for (int j = 0; j < A[0].size(); j++) {
				temp_value += A[i][j] * B[j][k];
			}
			result[i][k] = temp_value;
		}
	}
	return result;	
}

std::vector<double> CMatrix::Moore_Penrose(std::vector<double> _R, std::vector<std::vector<double>> J) {
	/// (J'*J)^-1 * J' * R
	double inv = 0.0;
	std::vector<double> R = _R;
	for (int i = 0; i < R.size(); i++)
		R[i] = -_R[i];
	std::vector<double> r_inv(R.size(), 0.0);
	std::vector<double> result(J.size(), 0.0);

	for (int i = 0; i < R.size(); i++)
		inv += R[i] * R[i];
	inv = 1 / inv;
	
	for (int i = 0; i < R.size(); i++)
		r_inv[i] = R[i] * inv;

	for (int i = 0; i < J.size(); i++)
			for (int k = 0; k < J.size(); k++)
				result[i] += R[k] * J[k][i];
			
	return result;
	}

std::vector<double> CMatrix::GaussSolver(std::vector<double> _R, std::vector<std::vector<double>> J) {
	/// para resolver o sistema linear, calculo a matriz inversa de J por eliminacao gaussiana, e multiplico por -R
	std::vector<double> result(J.size(), 0.0);
	std::vector<std::vector<double>> inv = CMatrix::inverse(J);

	for (int i = 0; i < J.size(); i++)
		for (int k = 0; k < J.size(); k++)
			result[i] += -_R[k] * inv[i][k];

	return result;
}

std::vector<std::vector<double>> CMatrix::inverse(std::vector<std::vector<double>> matrix) {
	int _size = matrix.size();
	std::vector<std::vector<double>> _inv(_size, std::vector<double>(_size));

	for (int i = 0; i < _size; i++) {
		_inv[i][i] = 1.0;
	}

	double pivot;
	// Gaussian elimination
	for (int p = 0; p < _size; p++) {
		// clean the row of pivot
		pivot = matrix[p][p];
		for (int i = 0; i < _size; i++) {
			matrix[p][i] = matrix[p][i] / pivot;
			_inv[p][i] = _inv[p][i] / pivot;
		}

		// clean the rows lower than pivot
		for (int i = 0; i < _size; i++) {
			pivot = matrix[i][p];
			if ((pivot < -1.0e-3 || pivot>1.0e-3) && i !=p) {
				for (int j = 0; j < _size; j++) {
					matrix[i][j] -= matrix[p][j] * pivot;
					_inv[i][j] -= _inv[p][j] * pivot;
				}
			}
		}
	}
	return _inv;
}

double CMatrix::determinanteSparse(std::vector<std::vector<double>> matriz){ //to find determinant {
	int n = matriz.size();
	double diagDireita = 1.0;
	double diagCentral = 1.0;
	double diagEsquerda = 1.0;

	for (int i = 0; i < n; i++) {
		diagCentral *= matriz[i][i];
		diagEsquerda *= (i == 0 ? 1.0 : matriz[i][i - 1]);
		diagDireita *= (i == n-1 ? 1.0 : matriz[i][i + 1]);
	}
	return diagDireita + diagCentral + diagEsquerda;
}

std::vector<double> CMatrix::LU_solver(std::vector<std::vector<double>> J, std::vector<double> R) {
	// decomposition of matrix
	size_t n = J.size();
	std::vector<std::vector<double>> lu(n);
	for (int i = 0; i < n; i++)
		lu[i].resize(n);

	double sum = 0;
	for (int i = 0; i < n; i++)	{
		for (int j = i; j < n; j++)	{
			sum = 0;
			for (int k = 0; k < i; k++)
				sum += lu[i][k] * lu[k][j];
			lu[i][j] = J[i][j] - sum;
		}
		for (int j = i + 1; j < n; j++) {
			sum = 0;
			for (int k = 0; k < i; k++)
				sum += lu[j][k] * lu[k][i];
			lu[j][i] = (1 / lu[i][i]) * (J[j][i] - sum);
		}
	}

	// lu = L+U-I
	// find solution of Ly = b
	std::vector<double> y(n, 0.0);
	for (int i = 0; i < n; i++)	{
		sum = 0;
		for (int k = 0; k < i; k++)
			sum += lu[i][k] * y[k];
		y[i] = -R[i] - sum; /// aqui tenho que subtrair R, pois a conta eh -R/J
	}
	// find solution of Ux = y
	std::vector<double> x(n, 0.0);

	for (int i = n - 1; i >= 0; i--) {
		sum = 0;
		for (int k = i + 1; k < n; k++)
			sum += lu[i][k] * x[k];
		x[i] = (1 / lu[i][i]) * (y[i] - sum);
	}
	return x;
}

double CMatrix::maxAbs(std::vector<double> R) {
	double max = 0.0;
	for (int i = 0; i < R.size(); i++)
		if (max < abs(R[i]))
			max = abs(R[i]);
	return max;
}

std::vector<std::vector<double>> CMatrix::rand(int n) {
	std::vector<std::vector<double>> _rand(n, std::vector<double>(n));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			_rand[i][j] = std::rand();
	return _rand;
}

void CMatrix::mostrarMatriz(std::vector<double> vector, std::string txt, bool stop) {
	std::cout << "\nMostrando o vetor " << txt << std::endl;
	for (unsigned int i = 0; i < vector.size(); i++)
		std::cout << std::setw(8) << std::setprecision(3) << vector[i] << " | ";

	if (stop) { char t[2]; std::cin.getline(t, 2); }
}

void CMatrix::mostrarMatriz(std::vector<std::vector<double>> matriz, std::string txt, bool stop) {
	std::cout << "\nMostrando a matriz " << txt << std::endl;
	for (unsigned int i = 0; i < matriz.size(); i++) {
		for (unsigned int k = 0; k < matriz[0].size(); k++) {
			std::cout << std::setw(8) << std::setprecision(3) << matriz[i][k] << " | ";
		}
		std::cout << std::endl;
	}
	if (stop) { char t[2]; std::cin.getline(t, 2); }
}