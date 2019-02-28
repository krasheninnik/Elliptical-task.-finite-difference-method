#include "pch.h"
#include "EllipticalTask.h"
#include "setting.h"
#include <iostream>

int main()
{
	EllipticalTask T;
	
#ifdef EXPLORE_CONVERGENCE
	const int amount_iters = 8;
	for(int i = 0; i < amount_iters; i++) {
		std::cout << "\riteration " << i << " of " << amount_iters << std::endl;
		T.init(i);

#ifdef CALCULATE
		T.solve(true);	// solve task
#else
		T.test();		// check the correctness of the decision on the given functions
#endif
	}

#else
	T.init(0);
#ifdef CALCULATE
	T.solve(true);	// solve task
#else
	T.test();		// check the correctness of the decision on the given functions
#endif
#endif
	
}
