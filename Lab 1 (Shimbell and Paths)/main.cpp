#include "Graph.h"

int main()
{
	Graph test;
	//test.Start();
	std::vector<std::vector<int>> explosion =
	{
		{0,0,12,0,51,0,0,0,0,0},
		{0,0,0,0,0,1,0,0,0,0},
		{-12,0,0,0,0,0,0,0,0,-62},
		{0,0,0,0,0,0,0,74,45,0},
		{0,0,0,0,0,0,0,-22,7,0},
		{0,0,0,0,0,0,-44,0,14,0},
		{0,0,0,0,0,0,0,52,-13,0},
		{0,0,0,0,0,0,0,0,52,-98},
		{0,0,67,0,0,0,0,0,0,-99},
		{0,0,62,0,0,0,0,98,99,0},
	};
	std::vector<int> path;
	test.bellmanFordForFlow(0, explosion, path);
	return 0;
}