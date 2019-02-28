#pragma once
#include <vector>
#include "grid.h"
#include "setting.h"
#include <functional>

class Matrix {
public:
	void init(Grid &G);
	void forming(Grid &G, double, double, double);

	int get_dim();

	std::vector<double> method_block_relaxation(std::vector<double>& x0, 
		std::vector<double>& f, std::vector<double>& temp, double w,
		std::vector<double>&, std::vector<double>&);

	void block_factorization(int, std::vector<double>&, std::vector<double>&);
	std::vector<double> solve_of_block_system(std::vector<double>, std::vector<double>, int,
		std::vector<double>&, std::vector<double>&);

	double calc_sum(int, std::vector<double>&);// sum of multiplicate elems of row matrix with corresponding vector's elems
	double calc_relative_discrepancy(std::vector<double> &, std::vector<double>&, std::vector<double>&);
	std::vector<double> multiplicate_with_vector(std::vector<double>&, std::vector<double>&);

	///////////
	std::vector<double> method_Jacobi(std::vector<double> &x1, std::vector<double> &x0,
		std::vector<double> &f, std::vector<double> &temp, double w);

private: 
	static const int diags = 5;

	std::vector<double> low_diag;

	std::vector<double> mid_diag1;
	std::vector<double> mid_diag2;
	std::vector<double> mid_diag3;

	std::vector<double> high_diag;
	
	int n, m, block_size, max_iter;
	double accuracy;
};

