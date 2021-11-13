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

void CWell::print() {
	std::cout << "\nObjeto Well" << std::endl;
	std::cout << "tp: ";
	for (int i = 0; i < tp.size(); i++) {
		std::cout << tp[i] << "  ";
	}

	std::cout << "\nqsc: ";
	for (int i = 0; i < qsc.size(); i++) {
		std::cout << qsc[i] << "  ";
	}
	std::cout << "\ndz | partial: \n";
	for (int i = 0; i < dz.size(); i++) {
		std::cout << dz[i] << " | " << partial[i] << std::endl;
	}
	std::cout << "rw: " << rw << std::endl;
}