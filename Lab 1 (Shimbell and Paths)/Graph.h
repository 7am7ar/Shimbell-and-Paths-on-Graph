#pragma once
#include <vector>

class Graph
{
public:
	Graph();
	void Start();
private:
	int dijkstra(int startVertex);
	void createGraph(int vertexQuantity);
	int bellmanFord(int startVertex);
	int floydWarshall(int startVertex);
	void shimbellMethod(int edgeQuantity, bool mode);
	void showMatrix(std::vector<std::vector<int>>& matrix);
	void showWeights(std::vector<std::vector<int>>& matrix);
	int dfs(int firstVertex, int secondVertex);

	int m_vertexQuantity;
	int m_mode;
	std::vector<std::vector<int>> m_matrix;
	std::vector<std::vector<int>> m_weightedMatrix;
	std::vector<int> m_outdegrees;
};