#include "CReservoir.hpp"

double CReservoir::calc_phi(double p) {
	return phi0 * (1.0 + cphi * (p - p0));
}

double CReservoir::calc_dphidp(double p) {
	return (phi0 * cphi);
}