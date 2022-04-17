#pragma once
#include <vector>

class Graph
{
public:
	Graph();
	void Start();
private:
	int dijkstra(int startVertex);
	int bellmanFord(int startVertex);
	void shimbellMethod(int edgeQuantity, bool mode);
	void showMatrix(std::vector<std::vector<int>>& matrix);
	void showWeights(std::vector<std::vector<int>>& matrix);
	int dfs(int firstVertex, int secondVertex);

	int m_vertexQuantity;
	std::vector<std::vector<int>> m_matrix;
	std::vector<std::vector<int>> m_weightedMatrix;
	std::vector<int> m_outdegrees;
};