#ifndef CWELL_HPP
#define CWELL_HPP

#include <vector>

class CWell {

private:
	std::vector<double> tp{0, 300};			/// tempos de mudanca na vazao na superficie [h]
	std::vector<double> qsc{0, -300 };		/// vazoes nos tempos de mudanca [m^3 std / dia]
	std::vector<double> partial{1,1}; /// porcentagem de abertura de cada delta z
	double rw{ 0.09486 };					/// raio do poco [m]

public:
	std::vector<double> get_tp() { return tp; }
	std::vector<double> get_qsc() { return qsc; }
	std::vector<double> get_partial() { return partial; }
	double get_rw() { return rw; }
	double get_qsc(double time);
	double get_partial(int i) { return partial[i]; }
};
#endif