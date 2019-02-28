
#include "pch.h"

#include "matrix.h"
#include "stdio.h"
#include <string>
#include <algorithm>
#include <math.h>
#include "conio.h"

double accuracy = 1e-15;

void Matrix::block_factorization(int iter, std::vector<double>& mid_factor_diag1, std::vector<double>& mid_factor_diag2) {	// LU decomposite
	if (mid_factor_diag2[iter * block_size] == 0) printf_s("error in LU decomposition. zero on diag.\n");

	for (int start_pos = iter * block_size; start_pos < (iter + 1) * block_size - 1; start_pos++) {
		//system("pause");
		mid_factor_diag1[start_pos] /= mid_factor_diag2[start_pos];

		if (abs(mid_factor_diag2[start_pos + 1] - mid_factor_diag1[start_pos] * mid_diag3[start_pos]) < mid_factor_diag2[start_pos + 1] * accuracy) {
			printf("error in LDU decomposition. zero on diag.\n");
		}
		else {
			mid_factor_diag2[start_pos + 1] -= mid_factor_diag1[start_pos] * mid_diag3[start_pos];
		}
	}
}

std::vector<double> Matrix::solve_of_block_system(std::vector<double> f, std::vector<double> y, int iter,
	std::vector<double>& mid_factor_diag1, std::vector<double>& mid_factor_diag2) {
	int i, matrix_index;

	// LUx = F
	// LY = F
	// Calculate Y:

	y[0] = f[0];

	for (i = 1, matrix_index = iter * block_size; i < block_size; i++, matrix_index++) {
		y[i] = f[i] - mid_factor_diag1[matrix_index] * y[i - 1];
	}

	// Ux = Y
	// Calculate x: (will storage x in vector f)

	f[block_size - 1] = y[block_size - 1] / mid_factor_diag2[(iter + 1) * block_size - 1];
	
	for (i = block_size - 2, matrix_index = (iter + 1) * block_size - 2; i >= 0; i--, matrix_index--) {
		f[i] = (y[i] - mid_diag3[matrix_index]*f[i+1]) / mid_factor_diag2[matrix_index];
	}

	return f;
}

std::vector<double> Matrix::method_block_relaxation(std::vector<double>& x0,
								std::vector<double>& f, std::vector<double>&temp, double w,
								std::vector<double>& mid_factor_diag1, std::vector<double>& mid_factor_diag2) {


	if (n % block_size != 0) {
		printf_s("Matrix don't divide into block's with size is %d\n", block_size);
		_getch();
		return x0;
	}

	if (block_size > m) {
		printf_s("block size too large: %d, max is %d\n", block_size, m);
		_getch();
		return x0;
	}

	std::vector<double> y(block_size),
		sum(block_size);

	int iter, i, k, line;
	double relative_discrepancy;

	// Factorization blocks of matrix
	mid_factor_diag1 = mid_diag1;
	mid_factor_diag2 = mid_diag2;

	for (i = 0; i < n / block_size; i++) {
		block_factorization(i, mid_factor_diag1, mid_factor_diag2);
	}
	

	for(iter = 1; iter < max_iter; iter++) {
		for (i = 0; i < n / block_size; i++) {	

			for (line = i * block_size, k = 0; k < block_size; line++, k++) {// Right parts.
				sum[k] = w*(f[line] - calc_sum(line, x0));			 // wRi(k)
			}

			// SOLVE OF SYSTEM AijYi(k) = wRi(k)
			y = solve_of_block_system(sum, y, i, mid_factor_diag1, mid_factor_diag2);
		
			// x1 = Y + x0
			for (line = i * block_size, k = 0; k < block_size; line++, k++) {
				x0[line] += y[k];
			}
		}

		// check discrepancy
		relative_discrepancy = calc_relative_discrepancy(x0, f, temp);

		int OUTPUT_COEF = max_iter / 100;
		if(iter % OUTPUT_COEF == 0) printf_s("\riteration: %d \t relative discrepancy: %e", iter, relative_discrepancy);

		if (relative_discrepancy < accuracy) {
			return x0;
		}
	}

	printf("i > max_iters!");
	return x0;		// return x? or error , amount iter grether then max_iter
}