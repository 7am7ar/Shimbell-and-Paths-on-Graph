#pragma once
#include <vector>

class Graph
{
public:
	Graph();
	void ShowMatrix();
private:
	int vertexQuantity;
	std::vector<std::vector<int>> matrix;
	std::vector<int> outdegrees;
};