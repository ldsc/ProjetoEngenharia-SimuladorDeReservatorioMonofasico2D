#include "CGas.hpp"

double CGas::calc_b(double p) {
	return (T0 * p) / (p0 * Z_KIAM(p) * T);
}

double CGas::calc_rho(double p) {
	return (p * Ma) / (Z_KIAM(p) * R * T);
}

double CGas::calc_dbdp(double p) {
	return T0 / (p0 * T * p * DELTA) * ((p * (1.0 + DELTA) / Z_KIAM(p * (1.0 + DELTA))) - (p / Z_KIAM(p)));
}

double CGas::calc_mu(double p) {
	return mu_LGE(0.001 * calc_rho(p));
}

double CGas::calc_dmudp(double p) {

	double pf = p * (1.0 + DELTA);
	double Zf = Z_KIAM(pf);

	double rhof = (pf * Ma) / (Zf * R * T);
	double muf = mu_LGE(0.001 * rhof);
	double mu = mu_LGE(0.001 * calc_rho(p));

	return (muf - mu) / (p * DELTA);
}

double CGas::Z_KIAM(double p) {
	double ppr = p / Ppc;
	double tpr = T / Tpc;

	double A1 = +0.3178420;
	double A2 = +0.3822160;
	double A3 = -7.7683540;
	double A4 = +14.290531;
	double A5 = +0.0000020;
	double A6 = -0.0046930;
	double A7 = +0.0962540;
	double A8 = +0.1667200;
	double A9 = +0.9669100;
	double A10 = +0.0630690;
	double A11 = -1.9668470;
	double A12 = +21.058100;
	double A13 = -27.024600;
	double A14 = +16.230000;
	double A15 = +207.78300;
	double A16 = -488.16100;
	double A17 = +176.29000;
	double A18 = +1.8845300;
	double A19 = +3.0592100;

	double t = 1 / tpr;
	double A = A1 * t * exp(A2 * pow(1.0 - t, 2)) * ppr;
	double B = A3 * t + A4 * pow(t, 2) + (A5 * pow(t, 6)) * pow(ppr, 6);
	double C = A9 + A8 * t * ppr + A7 * pow(t, 2) * pow(ppr, 2) + A6 * pow(t, 3) * pow(ppr, 3);
	double D = A10 * t * exp(A11 * pow(1.0 - t, 2));
	double E = A12 * t + A13 * pow(t, 2) + A14 * pow(t, 3);
	double F = A15 * t + A16 * pow(t, 2) + A17 * pow(t, 3);
	double G = A18 + A19 * t;
	double y = (D * ppr) / ((1 + pow(A, 2)) / C - (pow(A, 2)*B) / pow(C, 3));

	double Z = D * ppr*(1 + y + pow(y, 2) - pow(y, 3))
		/ ((D * ppr + E * pow(y, 2) - F * pow(y, G)) * pow(1 - y, 3));

	return Z;
}

double CGas::mu_LGE(double rho) {
	double K = ((9.379 + 0.01607 * Ma) * pow(1.8 * T, 1.5)) / (209.2 + 19.26 * Ma + 1.8 * T);
	double X = 3.448 + (986.4 / (1.8 * T)) + 0.01009 * Ma;
	double Y = 2.447 - 0.2224 * X;
	return (1.0e-4) * K * exp(X * pow(rho, Y));
}