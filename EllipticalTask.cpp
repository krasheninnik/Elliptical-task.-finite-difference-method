#include "pch.h"
#include "EllipticalTask.h"
#include <ostream>
#include <iomanip>

void EllipticalTask::formingRightPart() {
	for (int i = 0; i < G.nodes.size(); i++) {
		if (G.nodes[i].type != -1) {
			f[i] = funcs[G.nodes[i].type](G.nodes[i]);
		}
		else {
			f[i] = 0;
		}
	}
} 


void EllipticalTask::init(int k) {
	// read info about coefficients:
	std::fstream fin("equal_coefs.txt");
	fin >> lambda >> gamma >> beta;
	fin.close();

	setBoundaryCondition();

	// Memory allocation for:
	G.load(k);				// Grid
	M.init(G);				// Matrix

	uint32_t size = M.get_dim();
	
	f.resize(size);			// vector of right part
	x0.resize(size);		// vector of initial approximation
	u.resize(size);			// vector of result
	temp.resize(size);		// temp vector
	mfd1.resize(size-1);	// vector for mid_fatorized_diag1
	mfd2.resize(size);;		// vector for mid_fatorized_diag2
}


#include <iostream>
void EllipticalTask::solve(bool calculate) {
	M.forming(G, lambda, gamma, beta);			// forming Matrix from Grid and task

	if(calculate) formingRightPart();		    // forming Right part


	double w = 1;
	u = M.method_block_relaxation(x0, f, temp, w, mfd1, mfd2); // solve
}

void EllipticalTask::setBoundaryCondition() {
	/*
	funcs.resize(9);		// array of functions 

	funcs[0] = [](const Node &N) {return 0; };		// f: right part of LU = f(x,y)
	
	// dU/dn = - dU/dy:
	funcs[1] = [](const Node &N) {return 0; };		// boundary conditions:
	funcs[3] = [](const Node &N) {return 0; };
	funcs[7] = [](const Node &N) {return 0; };

	// dU/dn = dU/dy
	funcs[5] = [](const Node &N) {return 0; };

	// dU/dn = dU/dx
	funcs[2] = [](const Node &N) {return 0; };
	funcs[4] = [](const Node &N) {return 0; };
	
	// dU/dn = -dU/dx
	funcs[6] = [](const Node &N) {return 0; };
	funcs[8] = [](const Node &N) {return 0; };
	*/
}

//  -div(Lambda*gradU) + Gamma*U = f
void EllipticalTask::test() {
	//func uExactFunc = [](const Node &N) {return exp(N.x + N.y); };		// u = U(x,y)
	//func uExactFunc = [](const Node &N) {return exp(N.x); };		// u = U(x,y)

	func uExactFunc = [](const Node &N) {return pow(N.x,4) + pow(N.y,4); };

	std::vector<func> testFuncs(9);
//#define EXP

#ifdef BOUNDARY1
#ifdef EXP
	testFuncs[0] = [&](const Node &N) {return -lambda *  (2 *exp(N.x + N.y)) + gamma * (exp(N.x + N.y)); };		// f: f = -div(Lambda*gradU) + Gamma*U 

	// dU/dn = - dU/dy:
	testFuncs[1] = [&](const Node &N) {return gamma * (exp(N.x + N.y)); };		// boundary conditions:
	testFuncs[3] = [&](const Node &N) {return gamma * (exp(N.x + N.y)); };
	testFuncs[7] = [&](const Node &N) {return gamma * (exp(N.x + N.y)); };

	// dU/dn = dU/dy
	testFuncs[5] = [&](const Node &N) {return gamma * (exp(N.x + N.y)); };

	// dU/dn = dU/dx
	testFuncs[2] = [&](const Node &N) {return gamma * (exp(N.x + N.y)); };
	testFuncs[4] = [&](const Node &N) {return gamma * (exp(N.x + N.y)); };

	// dU/dn = -dU/dx
	testFuncs[6] = [&](const Node &N) {return gamma * (exp(N.x + N.y)); };
	testFuncs[8] = [&](const Node &N) {return gamma * (exp(N.x + N.y)); };
#else
	testFuncs[0] = [&](const Node &N) {return -lambda * (12 * pow(N.x,2) + 12 * pow(N.y, 2)) + gamma * (pow(N.x, 4) + pow(N.y, 4)); };

// dU/dn = - dU/dy:
	testFuncs[1] = [&](const Node &N) {return gamma * (pow(N.x, 4) + pow(N.y, 4)); };		// boundary conditions:
	testFuncs[3] = [&](const Node &N) {return gamma * (pow(N.x, 4) + pow(N.y, 4)); };
	testFuncs[7] = [&](const Node &N) {return gamma * (pow(N.x, 4) + pow(N.y, 4)); };
													  
	// dU/dn = dU/dy								 
	testFuncs[5] = [&](const Node &N) {return gamma * (pow(N.x, 4) + pow(N.y, 4)); };
													  
	testFuncs[2] = [&](const Node &N) {return gamma * (pow(N.x, 4) + pow(N.y, 4)); };
	testFuncs[4] = [&](const Node &N) {return gamma * (pow(N.x, 4) + pow(N.y, 4)); };
													  
	// dU/dn = -dU/dx								  
	testFuncs[6] = [&](const Node &N) {return gamma * (pow(N.x, 4) + pow(N.y, 4)); };
	testFuncs[8] = [&](const Node &N) {return gamma * (pow(N.x, 4) + pow(N.y, 4)); };
#endif
#else
#ifdef BOUNDARY2

	testFuncs[0] = [&](const Node &N) {return gamma * (pow(N.x,2)); };	// f: f = -div(Lambda*gradU) + Gamma*U 

	// dU/dn = - dU/dy:
	testFuncs[1] = [&](const Node &N) {return lambda * 0; };			// boundary conditions:
	testFuncs[3] = [&](const Node &N) {return lambda * 0; };
	testFuncs[7] = [&](const Node &N) {return lambda * 0; };

	// dU/dn = dU/dy
	//testFuncs[5] = [&](const Node &N) {return lambda * 0; };
	testFuncs[5] = [&](const Node &N) {return (pow(N.x, 2)); };


	// dU/dn = dU/dx
	testFuncs[2] = [&](const Node &N) {return lambda * 2 * N.x; };
	testFuncs[4] = [&](const Node &N) {return lambda * 2 * N.x; };

	// dU/dn = -dU/dx
	testFuncs[6] = [&](const Node &N) {return lambda * -1 * 2 * N.x; };
	testFuncs[8] = [&](const Node &N) {return lambda * -1 * 2 * N.x; };
	
#else
#ifdef BOUNDARY3	// lambda * (d/dn U

	testFuncs[0] = [&](const Node &N) {return -lambda * (12 * pow(N.x,2) + 12 * pow(N.y, 2)) + gamma * (pow(N.x, 4) + pow(N.y, 4)); };
	//testFuncs[0] = [&](const Node &N) {return -lambda*(6*N.x + 6*N.y) + gamma*(pow(N.x, 3) + pow(N.y, 3)); };			// f: f = -div(Lambda*gradU) + Gamma*U


	// dU/dn = - dU/dy:
	testFuncs[1] = [&](const Node &N) {return lambda * (-4 * pow(N.y, 3)) +  beta * (pow(N.x, 4) + pow(N.y, 4)); };		// boundary conditions:
	testFuncs[3] = [&](const Node &N) {return lambda * (-4 * pow(N.y, 3)) + beta * (pow(N.x, 4) + pow(N.y, 4)); };
	testFuncs[7] = [&](const Node &N) {return lambda * (-4 * pow(N.y, 3)) + beta * (pow(N.x, 4) + pow(N.y, 4)); };

	// dU/dn = dU/dy
	testFuncs[5] = [&](const Node &N) {return lambda * (4 * pow(N.y, 3)) + beta * (pow(N.x, 4) + pow(N.y, 4)); };

	// dU/dn = dU/dx
	testFuncs[2] = [&](const Node &N) {return lambda * (4 * pow(N.x, 3)) + beta * (pow(N.x, 4) + pow(N.y, 4)); };
	testFuncs[4] = [&](const Node &N) {return lambda * (4 * pow(N.x, 3)) + beta * (pow(N.x, 4) + pow(N.y, 4)); };

	// dU/dn = -dU/dx
	testFuncs[6] = [&](const Node &N) {return lambda * (-4 * pow(N.x, 3)) + beta * (pow(N.x, 4) + pow(N.y, 4)); };
	testFuncs[8] = [&](const Node &N) {return lambda * (-4 * pow(N.x, 3)) + beta * (pow(N.x, 4) + pow(N.y, 4)); };
#endif
#endif
#endif

	std::vector<double> uExact(M.get_dim());

	G.calculateValues(uExact, f, uExactFunc, testFuncs);
	
	solve(false);

	/* calcSum debug
	std::vector<double> x(M.get_dim(), 1);

	std::cout << std::endl;
	for (int i = 0; i < M.get_dim(); i++) {
		uExact[i] = M.calc_sum(i, x);
		std::cout << i << " : " <<uExact[i] << std::endl;
	}
	std::cout << std::endl;
	*/

	 // Debug result vizualization:
	 
#ifdef DISPLAY	
	std::cout << "\nReal:\n";
	std::cout.precision(7);
	int WIDTH = 14;

	for (uint32_t i = G.height - 1; i != UINT32_MAX; i--) {
		for (uint32_t j = 0, offset = i * G.width; j < G.width; j++) {
			std::cout << std::setw(WIDTH)  << uExact[offset + j];
		}
		std::cout << std::endl;
	}
	std::cout << "\nNumerical:\n";
	for (uint32_t i = G.height - 1; i != UINT32_MAX; i--) {
		for (uint32_t j = 0, offset = i * G.width; j < G.width; j++) {
			std::cout << std::setw(WIDTH) << u[offset + j];
		}
		std::cout << std::endl;
	}

	std::cout << "\nError:\n";
	for (uint32_t i = G.height - 1; i != UINT32_MAX; i--) {
		for (uint32_t j = 0, offset = i * G.width; j < G.width; j++) {
			std::cout << std::setw(WIDTH) << u[offset + j] - uExact[offset + j];
		}
		std::cout << std::endl;
	}
	
	/*
	for (uint32_t i = 0; i < G.nodes.size(); i++) {
		std::cout << G.nodes[i].x << "\t " << G.nodes[i].y << std::endl;
	}
	*/
#endif


	for (uint32_t i = 0; i < G.width; i++) std::cout << G.nodes[i].x << " ";
	std::cout << std::endl;
	for (uint32_t i = 0; i < G.height; i++) std::cout << G.nodes[i*G.width].y << " ";
	std::cout << std::endl;



	uint32_t maxNode;
	double maxRelevantError = 0;
	double averageRelevantError = 0;

	std::vector<double> absError(G.nodes.size());
	
	double sum = 0;
	for (uint32_t i = 0; i < G.nodes.size(); i++) sum += pow(u[i] - uExact[i], 2);
	sum = sqrt(sum);

	std::ofstream fout1("exploreGRID.txt", std::ios::app);
	fout1 << G.nodes[1].x - G.nodes[0].x << "\t" << sum << std::endl;
	fout1.close();




	std::ofstream fout("result.txt");
	fout.precision(16);
	//fout.width(20);
	//fout.fill(' ');
	fout << "node\tx\t   y\t\t\t\tnumerical\t\t\texact\t\t\terror" << std::endl;
	for (uint32_t i = 0; i < G.nodes.size(); i++) {
		fout << i << "\t" <<
			 std::setprecision(3) << std::scientific << G.nodes[i].x << "\t " <<
			std::scientific << G.nodes[i].y << "\t " <<
			std::setprecision(16) << std::scientific << u[i] << "\t" <<
			    std::scientific << uExact[i] << "\t" <<
				std::scientific << u[i] - uExact[i] << std::endl;
	}
}