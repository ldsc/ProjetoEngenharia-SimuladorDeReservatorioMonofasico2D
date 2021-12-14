#ifndef CFLUIDO_HPP
#define CFLUIDO_HPP

#include <string>

class CFluido {
public:
	virtual double calc_b(double p) { return 0.0; }
	virtual double calc_dbdp(double p) { return 0.0; }

	virtual double calc_mu(double p) { return 0.0; }
	virtual double calc_dmudp(double p) { return 0.0; }

	virtual std::string get_type() { return type; }
private:
	std::string type;

};

#endif