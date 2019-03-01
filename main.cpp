#include "pch.h"
#include "EllipticalTask.h"
#include "setting.h"
#include <iostream>

int main()
{
	EllipticalTask T;
	
#ifdef EXPLORE_CONVERGENCE
	const int amount_iters = 5;
	for(int i = 0; i < amount_iters; i++) {
		std::cout << "\riteration " << i << " of " << amount_iters << std::endl;
		T.init(i);
		T.test();		// check the correctness of the decision on the given functions
	}

#else
#ifdef CALCULATE
	T.init(0);
	T.solve(true);	// solve task
#else
	T.init(0);
	T.test();
#endif
#endif
	
}
