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
	bool dfsFordFulkerson(std::vector<std::vector<int>>& graph, std::vector<int>& currentPath, int startVertex, std::vector<bool>& isVisited);
	int fordFulkerson(int startVertex);
	void minCostFlow(int startVertex, int streamSize);
	bool bellmanFordForFlow(int startVertex, std::vector<std::vector<int>>& graph, std::vector<int>& path);
	bool findPath(std::vector<int>& path, int startVertex, std::vector<int>& marks, int finalVertex);
	bool findPathForFlow(std::vector<int>& path, int startVertex, std::vector<int>& marks, int finalVertex, std::vector<std::vector<int>>& graph, int lastVertex);
	bool isAchievable(int vertexOne, int vertexTwo, std::vector<std::vector<int>>& graph, int lastVertex);
	int findNumberOfSpanningTrees();
	int findDeterminant(std::vector<std::vector<int>>& matrix);
	int prim();
	int kruskal();
	void codePrufer();
	void decodePrufer();
	int m_vertexQuantity;
	int m_mode;
	std::vector<std::vector<int>> m_matrix;
	std::vector<std::vector<int>> m_weightedMatrix;
	std::vector<std::vector<int>> m_bandwidthMatrix;
	std::vector<std::vector<int>> m_minimumSpanningTree;
	std::vector<int> m_outdegrees;
	std::vector<std::pair<int, int>> m_pruferSpanningTree;
};