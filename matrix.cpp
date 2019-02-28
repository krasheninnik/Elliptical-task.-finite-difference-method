
#include "pch.h" 
#include "matrix.h"
#include "stdio.h"
#include <string>
#include <algorithm>
#include <math.h>
#include "conio.h"
#include <fstream>

double calc_norm(std::vector<double> &x) {
	double norm = 0;
	for (int i = 0; i < x.size(); i++) {
		norm += x[i] * x[i];
	}
	norm = sqrt(norm);
	return norm;
}

std::vector<double> load_vector(int size, std::string filename, std::vector<double> &vect) {
	FILE *in;
	fopen_s(&in, filename.c_str(), "r");
	for (int i = 0; i < size; i++) fscanf_s(in, "%lf ", &vect[i]);
	fclose(in);

	return vect;
}

void out_vector(std::vector<double> &x) {
	FILE *in;
	fopen_s(&in, "x.txt", "w");
	for (int i = 0; i < x.size(); i++) {
		fprintf_s(in, "lf", x[i]);
	}
	fclose(in);
}

void Matrix::init(Grid &G) {
	n = G.nodes.size();
	m = G.width;
	block_size = G.width;
	
	// read info about accuracy of decision
	std::fstream fin(R"(input\accuracy.txt)");
	fin >> max_iter >> accuracy;
	fin.close();

	// memory allocation for Matrix
	low_diag = std::vector<double>(n - G.width);
	mid_diag1 = std::vector<double>(n - 1);
	mid_diag2 = std::vector<double>(n);
	mid_diag3 = std::vector<double>(n - 1);
	high_diag = std::vector<double>(n - G.width);
}

void Matrix::forming(Grid &G, double lambda, double gamma, double beta){
	// calculate elems of matrix
	double hx1, hx2, hy1, hy2; // step for x and y

#ifdef CALCULATE	// calculate the task:

#else
#ifdef BOUNDARY1	// only first boundary conditionals: 
	for (int i = 0; i < G.nodes.size(); i++) {
		switch (G.nodes[i].type) {

		case -1:	//  fictitious node
			mid_diag2[i] = 1;
			break;

		case 0:
			hx1 = G.nodes[i].x - G.nodes[i - 1].x;	// hx(i-1)
			hx2 = G.nodes[i + 1].x - G.nodes[i].x;	// hx(i)
			hy1 = G.nodes[i].y - G.nodes[i - m].y;	// hy(i-1)
			hy2 = G.nodes[i + m].y - G.nodes[i].y;	// hy(i)

			low_diag[i - m] = -lambda * 2.0 / (hy1*(hy2 + hy1));	// miss the lamdba
			mid_diag1[i - 1] = -lambda * 2.0 / (hx1*(hx2 + hx1));
			mid_diag2[i] = -lambda * -2.0 * (1.0 / (hx1*hx2) + 1.0 / (hy1*hy2)) + gamma;
			mid_diag3[i] = -lambda * 2.0 / (hx2*(hx2 + hx1));
			high_diag[i] = -lambda * 2.0 / ((hy2*(hy2 + hy1)));
			break;

		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:
			mid_diag2[i] = 1;
			break;
		}
	}
#else
#ifdef BOUNDARY2	// only second boundary conditionals: 
	for (int i = 0; i < G.nodes.size(); i++) {
		switch (G.nodes[i].type) {
		case -1:	//  fictitious node
			mid_diag2[i] = 1;
			break;

		case 0:
			hx1 = G.nodes[i].x - G.nodes[i - 1].x;	// hx(i-1)
			hx2 = G.nodes[i + 1].x - G.nodes[i].x;	// hx(i)
			hy1 = G.nodes[i].y - G.nodes[i - m].y;	// hy(i-1)
			hy2 = G.nodes[i + m].y - G.nodes[i].y;	// hy(i)

			low_diag[i - m] = -lambda * 2.0 / (hy1*(hy2 + hy1)); // miss the lamdba
			mid_diag1[i - 1] = -lambda * 2.0 / (hx1*(hx2 + hx1));
			mid_diag2[i] = -lambda * -2.0*(1.0 / (hx1*hx2) + 1.0 / (hy1*hy2)) + gamma;
			mid_diag3[i] = -lambda * 2.0 / (hx2*(hx2 + hx1));
			high_diag[i] = -lambda * 2.0 / ((hy2*(hy2 + hy1)));
			break;


		case 1:	//mid_diag2[i] = 1; break;
		case 3: 	// dU/dn == - dU/dy
		case 7:
			hy2 = G.nodes[i + m].y - G.nodes[i].y;

			mid_diag2[i] =   lambda * 1.0/hy2;
			high_diag[i] = lambda *  -1.0/hy2;
			break;

		case 5:
			mid_diag2[i] = 1;
		
			/*
			hy1 = G.nodes[i].y - G.nodes[i - m].y;

			mid_diag2[i] = lambda * 1.0 / hy1;
			low_diag[i - m] = lambda * -1.0 / hy1;
			*/
			break;

		case 2:
		case 4:
			hx1 = G.nodes[i].x - G.nodes[i - 1].x;
			
			mid_diag2[i]   = lambda * 1.0 / hx1;
			mid_diag1[i-1] = lambda * -1.0 / hx1;
			break;

		case 6:	
		case 8:
			hx2 = G.nodes[i + 1].x - G.nodes[i].x;	

			mid_diag2[i] = lambda * 1.0 / hx2;
			mid_diag3[i] = lambda * -1.0 / hx2;
			break;
		}
	}
#else
#ifdef BOUNDARY3 // only third boundary conditionals: 
		for (int i = 0; i < G.nodes.size(); i++) {
			switch (G.nodes[i].type) {
			case -1:	//  fictitious node
				mid_diag2[i] = 1;
				break;

			case 0:
				hx1 = G.nodes[i].x - G.nodes[i - 1].x;	// hx(i-1)
				hx2 = G.nodes[i + 1].x - G.nodes[i].x;	// hx(i)
				hy1 = G.nodes[i].y - G.nodes[i - m].y;	// hy(i-1)
				hy2 = G.nodes[i + m].y - G.nodes[i].y;	// hy(i)

				low_diag[i - m] = -lambda * 2.0 / (hy1*(hy2 + hy1)); // miss the lamdba
				mid_diag1[i - 1] = -lambda * 2.0 / (hx1*(hx2 + hx1));
				mid_diag2[i] = -lambda * -2.0*(1.0 / (hx1*hx2) + 1.0 / (hy1*hy2)) + gamma;
				mid_diag3[i] = -lambda * 2.0 / (hx2*(hx2 + hx1));
				high_diag[i] = -lambda * 2.0 / ((hy2*(hy2 + hy1)));
				break;


			case 1:
			case 3: // dU/dn == - dU/dy
			case 7:
				hy2 = G.nodes[i + m].y - G.nodes[i].y;

				mid_diag2[i] = lambda * 1.0 / hy2 + beta ;
				high_diag[i] = lambda * -1.0 / hy2;
				break;

			case 5:

				hy1 = G.nodes[i].y - G.nodes[i - m].y;

				mid_diag2[i] = lambda * 1.0 / hy1 + beta;
				low_diag[i - m] = lambda * -1.0 / hy1;
				break;

			case 2:
			case 4:
				hx1 = G.nodes[i].x - G.nodes[i - 1].x;

				mid_diag2[i] = lambda * 1.0 / hx1 + beta;
				mid_diag1[i - 1] = lambda * -1.0 / hx1;
				break;

			case 6:
			case 8:
				hx2 = G.nodes[i + 1].x - G.nodes[i].x;

				mid_diag2[i] = lambda * 1.0 / hx2 + beta;
				mid_diag3[i] = lambda * -1.0 / hx2;
				break;
			}
		}
#endif
#endif
#endif
#endif
}

int Matrix::get_dim() { return n; }

// should work. need test calc_sum
double Matrix::calc_sum(int row, std::vector<double> &x) {
	// m - wight of grid
	int ldi = row - m;  // first low_diags_index
	int hdi = row + m;	// first high diags_index
	int i;
	double sum = 0;

	// Throw low diags
	if (ldi >= 0) { // use only 3-th low diags.
		sum += low_diag[ldi] * x[ldi];
	}


	// Throw main diags
	if (row > 0) {
		if (row < n - 1) { // use all main diags
			sum += mid_diag1[row - 1] * x[row - 1];		
			sum += mid_diag2[row] * x[row];				
			sum += mid_diag3[row] * x[row + 1];			
		}

		else {		// use 1,2 main diags
			sum += mid_diag1[row - 1] * x[row - 1];		
			sum += mid_diag2[row] * x[row];		
		}
	}
	else {			// use 2,3 main diags
		sum += mid_diag2[row] * x[row];	
		sum += mid_diag3[row] * x[row + 1];	
	}

	// Throw high diags 
	if (n - hdi > 0) {
		sum += high_diag[row] * x[hdi];
	}
	
	return sum;
}

std::vector<double> Matrix::multiplicate_with_vector(std::vector<double> &x, std::vector<double> &f) {
	for (int i = 0; i < n; i++) f[i] = calc_sum(i, x);
	return f;
}

double Matrix::calc_relative_discrepancy(std::vector<double> &x, std::vector<double> &f, std::vector<double> &ax) {
	ax = multiplicate_with_vector(x, ax);	    	      // Calculate Ax.
	for (int i = 0; i < ax.size(); i++) ax[i] -= f[i];    // Calculate ( Ax - f )
	double relative_discr = calc_norm(ax) / calc_norm(f); // Calculate discrepancy

	return relative_discr;
}

///////////
std::vector<double> Matrix::method_Jacobi(std::vector<double> &x1, std::vector<double> &x0,
	std::vector<double> &f, std::vector<double> &temp, double w) {
	printf_s("Method Jacobi: ");
	double relative_discrepancy;

	for (int i = 1; i <= max_iter; i++) {

		for (int j = 0; j < n; j++) {		// Calculate Xi of next iteration
			x1[j] = x0[j] + (w / mid_diag2[j]) * (f[j] - calc_sum(j, x0));
		}

		// check discrepancy
		relative_discrepancy = calc_relative_discrepancy(x1, f, temp);
		printf_s("\riteration: %d \t relative discrepancy: %e \n", i, relative_discrepancy);

		if (relative_discrepancy < accuracy) {
			return x1;
		}
		std::swap(x1, x0);		// preparate next iteration			
	}

	printf("i > max_iters!");

	return x0;		// return x? or error , amount iter grether then max_iter

}