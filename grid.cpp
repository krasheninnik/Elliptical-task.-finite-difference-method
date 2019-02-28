#include "pch.h"
#include "grid.h"
#include "setting.h"
#include <fstream>
#include <assert.h>

void Grid::load(int k) {
	using namespace std;
	fstream fin("grid_in.txt");

	uint32_t firstColumn;
	uint32_t secondColumn;
	uint32_t firstLine;

	// read x0 and y0
	double x0, y0;
	fin >> x0 >> y0;

	// input struct of grid
	fin >> firstColumn >> secondColumn >> width;
	fin >> firstLine >> height;

	// read x and y steps.
	double stepX, stepY;  
	fin >> stepX >> stepY;

	// consider grid's unevenness
	uint32_t amount;
	double coef;

#ifdef EXPLORE_CONVERGENCE
	uint32_t iter = pow(2, k); // 8
	firstColumn *= iter;
	secondColumn *= iter;
	width *= iter;
	firstLine *= iter;
	height *= iter;
	//stepX /= double(iter);
	//stepY /= double(iter);

#endif



	std::vector<std::pair<uint32_t, double>> unevennessX;
	for (int i = 0, count = 0; count != width; i++) {
		fin >> amount >> coef;
#ifdef EXPLORE_CONVERGENCE
		amount *= iter;	
		if (coef != 1) coef = pow(coef, 1.0 / iter);	// for unevenness grid need change coef too.
#endif
		unevennessX.push_back(std::pair<int, double>(amount, coef));
		count += amount;
	}

	std::vector<std::pair<uint32_t, double>> unevennessY;
	for (int i = 0, count = 0; count != width; i++) {
		fin >> amount >> coef;
#ifdef EXPLORE_CONVERGENCE
		amount *= iter;
		if (coef != 1) coef = pow(coef, 1.0 / iter);	// for unevenness grid need change coef too.
		
#endif
		unevennessY.push_back(std::pair<int, double>(amount, coef));
		count += amount;
	}

	///////////////////////////
	if(k != 0) {
		stepX /= 1 + unevennessX.back().second;
		stepY /= 1 + unevennessY.back().second;
	}
	////////////////////////


	// formation grig:
	// memory allocation for grid
	nodes = std::vector<Node>(width*height);		// can throw exception...

	double dx, dy;
	std::vector<double> dxs(width);
	
	size_t i = 0;

	for (int j = 0; j < unevennessX.size(); j++) {	// fill dxs array
		coef = unevennessX[j].second;				// coefficient of unevenness
		for (int k = 0; k < unevennessX[j].first; k++) {
			dxs[i++] = stepX;
			stepX *= coef;
		}
	}
		
	assert(i == dxs.size());

	double currentY = y0;


	for (int i = 0; i < unevennessY.size(); i++) {	// fill dxs array
		coef = unevennessX[i].second;				// coefficient of unevenness
		for (int k = 0; k < unevennessX[i].first; k++) {
			int offset = k * width;

			nodes[offset].x = x0;
			nodes[offset].y = currentY;

			for (int j = 1; j < width; j++) {
				nodes[offset + j].x = nodes[offset + j - 1].x + dxs[j - 1];
				nodes[offset + j].y = currentY;
			}

			currentY += stepY;
			stepY *= coef;
		}
	}

	// Calculate types of nodes
	/*   struct of grid:

		 55555555555555			//	0: internal node
		 60000000000004			//  N: node of boundary N
		 60000000000004			//  *: fictitious node
		 77777000033333
		 *****80002****
		 *****80002****
		 *****11111****
	*/
	 

	uint32_t offset;

	// mark fictitious nodes ++ 
	for (uint32_t i = 0; i < firstLine - 1; i++) {
		offset = i * width;

		for (uint32_t j = 0; j < firstColumn - 1; j++) {
			nodes[offset + j].type = -1;
		}

		for (uint32_t j = secondColumn; j < width; j++) {
			nodes[offset + j].type = -1;
		}
	}

	// mark internal nodes
	for (uint32_t i = 1; i < firstLine; i++) {
		offset = i * width;

		for (uint32_t j = firstColumn; j < secondColumn - 1; j++) {
			nodes[offset + j].type = 0;
		}
	}

	for (uint32_t i = firstLine; i < height - 1; i++) {
		offset = i * width;

		for (uint32_t j = 1; j < width - 1; j++) {
			nodes[offset + j].type = 0;
		}
	}

	// mark nodes of 1 boundary
	for (uint32_t j = firstColumn - 1; j < secondColumn; j++) {
		nodes[j].type = 1;
	}

	// mark nodes of 2 and 8 boundary
	for (uint32_t i = 1; i < firstLine - 1; i++) {
		offset = i * width;
		nodes[offset + firstColumn - 1].type = 8;
		nodes[offset + secondColumn - 1].type = 2;
	}

	// mark nodes of 7 and 3 boundary
	offset = (firstLine - 1) * width;

	for (uint32_t j = 0; j < firstColumn; j++) {
		nodes[offset + j].type = 7;
	}

	for (uint32_t j = secondColumn - 1; j < width; j++) {
		nodes[offset + j].type = 3;
	}

	// mark nodes of 6 and 4 boundary
	for (uint32_t i = firstLine; i < height - 1; i++) {
		offset = i * width;
		nodes[offset].type = 6;
		nodes[offset + width - 1].type = 4;
	}

	// mark nodes of 5 boundary
	offset = (height - 1) * width;

	for (uint32_t j = 0; j < width; j++) {
		nodes[offset + j].type = 5;
	}
}
/*
	Calculate exact values of the function  U and function of right part f on the grid.
	Uses for testig.
*/
void Grid::calculateValues(std::vector<double> &uExact, std::vector<double> &f,
							const  std::function<double(const Node&)> uExactFunc,
	const std::vector<std::function<double(const Node&)>> testFuncs) {

	for (int i = 0; i < nodes.size(); i++) {
		if (nodes[i].type != -1) {	// processing internal and boundary nodes
			uExact[i] = uExactFunc(nodes[i]);
			f[i] = testFuncs[nodes[i].type](nodes[i]);
		}
		else {					// processing fictious nodes
			uExact[i] = 0;
			f[i] = 0;
		}
		
	}
}