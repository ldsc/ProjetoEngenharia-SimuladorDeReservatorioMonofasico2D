#include "CWell.hpp"

double CWell::get_qsc(double t) {
	bool end = true;
	double q;
	for (int i = 0; i < tp.size(); i++) {
		if (t <= tp[i]) { // esse 1 eh so para gerar uma margem, eh muito literal o igual
			q = qsc[i];
			end = false;
		}
	}
	return end ? qsc[qsc.size() - 1] : q;
}