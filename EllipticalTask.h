#pragma once
#include "matrix.h"
#include "grid.h"
#include <fstream>
#include "setting.h"


typedef std::function<double(const Node&)> func;

class EllipticalTask {				// div (-lambda * grad U ) + gamma*U = f	
public:
	void solve(bool);
	void init(int);

	void test();

private:
	void formingRightPart();
	void setBoundaryCondition();	

	double lambda;					// coefficients for equal:	
	double gamma;
	double beta;


	Grid G;							// grid
	Matrix M;						// matrix
	std::vector<func> funcs;		// boundary conditional

	std::vector<double> f;			// vector of right part
	std::vector<double> x0;			// vector of initial approximation
	std::vector<double> u;			// vector of result
	std::vector<double> temp;		// temp vector
	std::vector<double> mfd1;		// vector for mid_fatorized_diag1
	std::vector<double> mfd2;		// vector for mid_fatorized_diag2
};