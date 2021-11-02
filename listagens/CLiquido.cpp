#include "CLiquido.hpp"

double CLiquido::calc_b(double p) {	return b0 * (1.0 + cf * (p - p0)); }

double CLiquido::calc_dbdp(double p) { return b0 * cf; }

double CLiquido::calc_mu(double p) { return mu * (1.0 + cmu * (p - p0)); }

double CLiquido::calc_dmudp(double p) {	return mu * cmu; }