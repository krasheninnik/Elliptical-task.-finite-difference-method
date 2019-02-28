#pragma once
#include "pch.h"
#include <functional>
#include <vector>

//typedef std::function<double(const Node&)> func;

struct Node {
	double x;
	double y;

	int type = -99;	/* type of Node: -1: fictitious
				//				  0: internal
				//				  n: number of boundary
				//				-99: undefined
				*/		
};

struct Grid {			// "T"-like grid
	void load(int);
	void calculateValues(std::vector<double>&, std::vector<double>&,
		const std::function<double(const Node&)>,
		const std::vector<std::function<double(const Node&)>>);

	uint32_t width;				
	uint32_t height;

	std::vector<Node> nodes;
};