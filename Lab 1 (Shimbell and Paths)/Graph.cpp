#include "Graph.h"
#include "Constants.h"
#include "GlobalFunctions.h"

#include <iostream>
#include <algorithm>
#include <random>
#include <ctime>
#include <set>
#include <queue>
#include <deque>
#include <list>

#define stop __asm nop

Graph::Graph()
{
	std::string input;
	bool isAssigned = false;
	while (!isAssigned)
	{
		std::cout << "Enter the number of vertexes:\n";
		std::cin >> input;
		if (IsOnlyDigits(input) && std::stoi(input) > 1) isAssigned = true;
		else std::cout << "Number is incorrect.\n";
	}
	createGraph(std::stoi(input));
}

void Graph::createGraph(int vertexQuantity)
{
	// Delete the previous graph
	m_matrix.clear();
	m_weightedMatrix.clear();
	m_bandwidthMatrix.clear();
	m_outdegrees.clear();
	m_vertexQuantity = vertexQuantity;
	m_minimumSpanningTree.clear();
	m_pruferSpanningTree.clear();

	//Calculate value of p0
	double p0 = 1;
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		p0 *= r + (i * c);
		p0 /= b + r + (i * c);
	}

	//Calculating m_outdegrees
	std::mt19937 mersenne(static_cast<unsigned int>(time(0)));
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		int x = 0;
		double p = p0;
		double random = static_cast<double>(mersenne() % 10001) / 10000;
		while (true)
		{
			random -= p;
			if (random < 0)
			{
				if (x == 0) x = 1;
				m_outdegrees.push_back(x);
				break;
			}
			else
			{
				if (x < m_vertexQuantity - 1)
				{
					x++;
					p = p * (b + c * (x - 1)) * (m_vertexQuantity - 1 - x);
					p = p / ((random + (m_vertexQuantity - x) * c) * x);
				}
				else
				{
					m_outdegrees.push_back(m_vertexQuantity);
					break;
				}
			}
		}
	}

	//Calculate actual m_matrix
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		std::vector<int> temp;
		//Fill with zeros values up to diagonal
		for (int j = 0; j < i + 1; j++) temp.push_back(0);

		//If Outdegree is more than empty space, fill with ones
		if (m_outdegrees.at(i) >= m_vertexQuantity - i - 1)
		{
			for (int j = i + 1; j < m_vertexQuantity; j++) temp.push_back(1);
		}
		else
		{
			std::set<int> pos;
			//Create set with positions of ones
			while (pos.size() != m_outdegrees.at(i))
			{
				pos.insert((mersenne() % (m_vertexQuantity - i - 1)) + i + 1);
			}
			//Fill empty space with appropriate values
			for (int j = i + 1; j < m_vertexQuantity; j++)
			{
				if (pos.find(j) != pos.end()) temp.push_back(1);
				else temp.push_back(0);
			}
		}
		m_matrix.push_back(temp);
	}
	m_matrix.at(m_vertexQuantity - 2).at(m_vertexQuantity - 1) = 1;

	//Calculating weights for edges
	std::string input;
	bool isAssigned = false;
	std::cout << "0. Only positive weights.\n";
	std::cout << "1. Only negative weights.\n";
	std::cout << "2. Mixed weights.\n";

	while (!isAssigned)
	{
		std::cout << "Enter the option number:\n";
		std::cin >> input;
		if (IsOnlyDigits(input) && std::stoi(input) >= 0 && std::stoi(input) <= 2) isAssigned = true;
		else std::cout << "Number is incorrect.\n";
	}
	m_mode = std::stoi(input);
	m_weightedMatrix = m_matrix;
	for (int i = 0; i < m_weightedMatrix.size(); i++)
	{
		for (int j = 0; j < m_weightedMatrix.at(i).size(); j++)
		{
			if (m_weightedMatrix.at(i).at(j) == 1) m_weightedMatrix.at(i).at(j) = (mersenne() % 100) + 1;
			if (m_mode == 1) m_weightedMatrix.at(i).at(j) *= -1;
			if (m_mode == 2) m_weightedMatrix.at(i).at(j) *= std::pow(-1, mersenne() % 2);
		}
	}

	//Calculate bandwidth
	m_bandwidthMatrix = m_matrix;
	for (int i = 0; i < m_bandwidthMatrix.size(); i++)
	{
		for (int j = 0; j < m_bandwidthMatrix.at(i).size(); j++)
		{
			if (m_bandwidthMatrix.at(i).at(j) == 1) m_bandwidthMatrix.at(i).at(j) = (mersenne() % 1000) + 1;
		}
	}
}


void Graph::showMatrix(std::vector<std::vector<int>>& matrix)
{
	std::cout << "\n\t";
	for (int i = 0; i < matrix.size(); i++)
	{
		std::cout << '(' << i << ")\t";
	}
	std::cout << '\n';

	for (int i = 0; i < matrix.size(); i++)
	{
		std::cout << '(' << i << ")\t";
		for (int j = 0; j < matrix.at(i).size(); j++)
		{
			std::cout << matrix.at(i).at(j) << '\t';
		}
		std::cout << '\n';
	}
}

void Graph::showWeights(std::vector<std::vector<int>>& matrix)
{
	bool isEmpty = true;
	std::cout << '\n';
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix.at(i).size(); j++)
		{
			if (matrix.at(i).at(j))
			{
				isEmpty = false;
				std::cout << i << "-->" << j << " Weight: " << matrix.at(i).at(j) << "\n";
			}
		}
	}
	if (isEmpty) std::cout << "There are no paths.";
}

void Graph::shimbellMethod(int edgeQuantity, bool mode)
{
	auto shimbellMatrix = m_weightedMatrix;
	auto tempMatrix = shimbellMatrix;

	for (int i = 0; i < edgeQuantity - 1; i++)
	{
		for (int j = 0; j < m_vertexQuantity; j++)
		{
			for (int k = 0; k < m_vertexQuantity; k++)
			{
				int currentValue = 0;
				for (int l = 0; l < m_vertexQuantity; l++)
				{
					if (shimbellMatrix.at(j).at(l) != 0 && m_weightedMatrix.at(l).at(k) != 0)
					{
						if (currentValue == 0) currentValue = shimbellMatrix.at(j).at(l) + m_weightedMatrix.at(l).at(k);
						else
						{
							if (mode) currentValue = std::min(currentValue, shimbellMatrix.at(j).at(l) + m_weightedMatrix.at(l).at(k));
							else currentValue = std::max(currentValue, shimbellMatrix.at(j).at(l) + m_weightedMatrix.at(l).at(k));
						}
					}
				}
				tempMatrix.at(j).at(k) = currentValue;
			}
		}
		shimbellMatrix = tempMatrix;
	}

	std::cout << "Shimbell method results:\n";
	showWeights(shimbellMatrix);
}

int Graph::dfs(int firstVertex, int secondVertex)
{
	int count = 0;
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		if (m_matrix.at(firstVertex).at(i) != 0)
		{
			if (i != secondVertex) count += dfs(i, secondVertex);
			else count++;
		}
	}
	return count;
}

int Graph::dijkstra(int startVertex)
{
	// Create start position
	std::vector<int> marks(m_vertexQuantity, INT_MAX);
	marks[startVertex] = 0;
	std::vector<int> isVisited(m_vertexQuantity, false);
	std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> queue;
	queue.push(std::make_pair(0, startVertex));
	int iterationCounter = 0;

	// Algorithm
	while (!queue.empty())
	{
		auto currentVertex = queue.top();
		queue.pop();
		isVisited[currentVertex.second] = true;
		for (int i = 0; i < m_vertexQuantity; i++)
		{
			iterationCounter++;
			if (m_weightedMatrix[currentVertex.second][i] != 0 && !(isVisited[i]))
			{
				if (marks[i] > marks[currentVertex.second] + m_weightedMatrix[currentVertex.second][i])
				{
					marks[i] = marks[currentVertex.second] + m_weightedMatrix[currentVertex.second][i];
					queue.push(std::make_pair(marks[i], i));
				}
			}
		}
	}

	//Output paths
	std::cout << "Dijkstra results:\n";
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		std::vector<int> path;
		if (marks[i] == INT_MAX) std::cout << startVertex << "-->" << i << " Shortest path length: " << "infinity" << '\n';
		else
		{
			findPath(path, startVertex, marks, i);
			std::cout << startVertex;
			while (!path.empty())
			{
				std::cout << "-->" << *(path.end() - 1);
				path.pop_back();
			}
			std::cout << "-->" << i << " Shortest path length: " << marks[i] << '\n';
		}
	}
	return iterationCounter;
}

bool Graph::findPath(std::vector<int>& path, int startVertex, std::vector<int>& marks, int finalVertex)
{
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		if (marks[i] + m_weightedMatrix[i][finalVertex] == marks[finalVertex] && m_weightedMatrix[i][finalVertex] != 0)
		{
			if (startVertex == i) return true;
			path.push_back(i);
			if (findPath(path, startVertex, marks, i))
			{
				return true;
			}
			else path.pop_back();
		}
	}
	return false;
}

int Graph::bellmanFord(int startVertex)
{
	//Create start position
	std::vector<int> marks(m_vertexQuantity, INT_MAX);
	std::vector<bool> isAppearedInQueue(m_vertexQuantity, false);
	std::vector<bool> isInQueue(m_vertexQuantity, false);
	marks[startVertex] = 0;
	isAppearedInQueue[startVertex] = true;
	std::deque<std::pair<int, int>> queue;
	queue.push_back(std::make_pair(0, startVertex));
	int iterationCounter = 0;

	//Algorithm
	while (!queue.empty())
	{
		auto currentVertex = queue[0];
		queue.pop_front();
		isInQueue[currentVertex.second] = false;

		for (int i = 0; i < m_vertexQuantity; i++)
		{
			iterationCounter++;
			if (m_weightedMatrix[currentVertex.second][i] != 0)
			{
				if (marks[i] > marks[currentVertex.second] + m_weightedMatrix[currentVertex.second][i])
				{
					marks[i] = marks[currentVertex.second] + m_weightedMatrix[currentVertex.second][i];
					if (isAppearedInQueue[i])
					{
						if (isInQueue[i])
						{
							for (auto j = queue.begin(); j != queue.end(); j++)
							{
								if ((*j).second == i)
								{
									queue.erase(j);
									break;
								}
							}
						}
						queue.push_front(std::make_pair(marks[i], i));
					}
					else
					{
						queue.push_back(std::make_pair(marks[i], i));
					}
					isAppearedInQueue[i] = true;
					isInQueue[i] = true;
				}
			}
		}
	}

	//Output marks
	std::cout << "Bellman-Ford results:\n";
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		std::vector<int> path;
		if (marks[i] == INT_MAX) std::cout << startVertex << "-->" << i << " Shortest path length: " << "infinity" << '\n';
		else
		{
			findPath(path, startVertex, marks, i);
			std::cout << startVertex;
			while (!path.empty())
			{
				std::cout << "-->" << *(path.end() - 1);
				path.pop_back();
			}
			std::cout << "-->" << i << " Shortest path length: " << marks[i] << '\n';
		}
	}
	return iterationCounter;
}

int Graph::floydWarshall(int startVertex)
{
	int iterationCounter = 0;
	std::vector<std::vector<int>> distance;
	distance = m_weightedMatrix;
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		for (int j = 0; j < m_vertexQuantity; j++)
		{
			if (j == i) distance[i][j];
			else if (distance[i][j] == 0) distance[i][j] = INT_MAX;
		}
	}

	for (int k = 0; k < m_vertexQuantity; k++)
	{
		for (int i = 0; i < m_vertexQuantity; i++)
		{
			for (int j = 0; j < m_vertexQuantity; j++)
			{
				iterationCounter++;
				if (distance[i][k] != INT_MAX && distance[k][j] != INT_MAX)
				{
					if (distance[i][j] > distance[i][k] + distance[k][j])
					{
						distance[i][j] = distance[i][k] + distance[k][j];
					}
				}
			}
		}
	}

	//Output distance
	std::cout << "Floyd-Warshall results:\n";
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		for (int j = 0; j < m_vertexQuantity; j++)
		{
			if (distance[i][j] == INT_MAX) std::cout << "inf";
			else std::cout << distance[i][j];
			std::cout << '\t';
		}
		std::cout << '\n';
	}
	return iterationCounter;
}

bool Graph::dfsFordFulkerson(std::vector<std::vector<int>>& graph, std::vector<int>& currentPath, int startVertex, std::vector<bool>& isVisited)
{
	bool result = false;

	for (int i = 0; i < m_vertexQuantity; i++)
	{
		if (graph[startVertex][i] > 0 && isVisited[i] == false)
		{
			currentPath.push_back(i);
			isVisited[startVertex] = true;
			if (i == m_vertexQuantity - 1) return true;
			result = dfsFordFulkerson(graph, currentPath, i, isVisited);
			if (result == false)
			{
				currentPath.pop_back();
				isVisited[startVertex] = false;
			}
			else return true;
		}
	}
	return false;
}

int Graph::fordFulkerson(int startVertex)
{
	std::cout << "Bandwidth matrix:";
	showMatrix(m_bandwidthMatrix);
	std::cout << '\n';

	int streamValue = 0;
	std::vector<std::vector<int>> streamMatrix(m_vertexQuantity, std::vector<int>(m_vertexQuantity, 0));
	std::vector<int> currentPath;
	std::vector<bool> isVisited(m_vertexQuantity, false);
	currentPath.push_back(startVertex);
	auto bandwidthMatrix = m_bandwidthMatrix;

	while (dfsFordFulkerson(bandwidthMatrix, currentPath, startVertex, isVisited))
	{
		int minWeight = INT_MAX;
		
		//Find maximum flow in path
		for (int i = 0; i < currentPath.size() - 1; i++)
		{
			if (bandwidthMatrix[currentPath[i]][currentPath[i + 1]] < minWeight)
				minWeight = bandwidthMatrix[currentPath[i]][currentPath[i + 1]];
		}
		streamValue += minWeight;
		//Change matrixes
		for (int i = 0; i < currentPath.size() - 1; i++)
		{
			bandwidthMatrix[currentPath[i]][currentPath[i + 1]] -= minWeight;
			bandwidthMatrix[currentPath[i + 1]][currentPath[i]] += minWeight;

			if (currentPath[i] < currentPath[i + 1]) 
				streamMatrix[currentPath[i]][currentPath[i + 1]] += minWeight;
			if ((streamMatrix[currentPath[i + 1]][currentPath[i]] - minWeight) >= 0) 
				streamMatrix[currentPath[i + 1]][currentPath[i]] -= minWeight;
		}

		//Find next path
		currentPath.clear();
		currentPath.push_back(startVertex);
		for (int i = 0; i < m_vertexQuantity; i++) isVisited[i] = false;
	}

	//Output stream values
	int count = 0;

	for (int i = 0; i < m_vertexQuantity; i++)
	{
		count += streamMatrix[i][m_vertexQuantity - 1];
	}
	std::cout << "Ford-Fulkerson results: ";
	std::cout << "\nMax flow: " << count << '\n';
	showMatrix(streamMatrix);
	std::cout << '\n';
	return streamValue;
}

bool Graph::findPathForFlow(std::vector<int>& path, int startVertex, std::vector<int>& marks, int finalVertex, std::vector<std::vector<int>>& graph, int lastVertex)
{
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		if (i != lastVertex)
		{
			if (marks[i] + graph[i][finalVertex] == marks[finalVertex] && graph[i][finalVertex] != 0)
			{
				if (startVertex == i) return true;
				path.push_back(i);
				if (findPathForFlow(path, startVertex, marks, i, graph, finalVertex))
				{
					return true;
				}
				else path.pop_back();
			}
		}
	}
	return false;
}

bool Graph::bellmanFordForFlow(int startVertex, std::vector<std::vector<int>>& graph, std::vector<int>& path)
{
	//Create start position
	std::vector<int> marks(m_vertexQuantity, INT_MAX);
	std::vector<bool> isAppearedInQueue(m_vertexQuantity, false);
	std::vector<bool> isInQueue(m_vertexQuantity, false);
	marks[startVertex] = 0;
	isAppearedInQueue[startVertex] = true;
	std::deque<std::pair<int, int>> queue;
	queue.push_back(std::make_pair(0, startVertex));

	//Algorithm
	while (!queue.empty())
	{
		auto currentVertex = queue[0];
		queue.pop_front();
		isInQueue[currentVertex.second] = false;

		for (int i = 0; i < m_vertexQuantity; i++)
		{
			if (i != currentVertex.first)
			{
				if (graph[currentVertex.second][i] != 0)
				{
					if (marks[i] > marks[currentVertex.second] + graph[currentVertex.second][i])
					{
						marks[i] = marks[currentVertex.second] + graph[currentVertex.second][i];
						if (i != m_vertexQuantity - 1)
						{
							if (isAppearedInQueue[i])
							{
								if (isInQueue[i])
								{
									for (auto j = queue.begin(); j != queue.end(); j++)
									{
										if ((*j).second == i)
										{
											queue.erase(j);
											break;
										}
									}
								}
								queue.push_front(std::make_pair(currentVertex.second, i));
							}
							else
							{
								queue.push_back(std::make_pair(currentVertex.second, i));
							}
							isAppearedInQueue[i] = true;
							isInQueue[i] = true;
						}
					}
				}
			}
		}
	}

	//Find the path to last vertex
	if (marks[m_vertexQuantity - 1] == INT_MAX)
	{
		return false;
	}
	else findPathForFlow(path, startVertex, marks, m_vertexQuantity - 1, graph, -1);
	return true;
}

void Graph::minCostFlow(int startVertex, int streamSize)
{
	int streamValue = 0;
	bool isOver = false;
	std::vector<std::vector<int>> streamMatrix(m_vertexQuantity, std::vector<int>(m_vertexQuantity, 0));
	std::vector<int> currentPath;
	currentPath.push_back(startVertex);
	auto bandwidthMatrix = m_bandwidthMatrix;
	auto weightedMatrix = m_weightedMatrix;
	auto modifiedWeightedMatrix = m_weightedMatrix;
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		for (int j = 0; j < i; j++)
		{
			modifiedWeightedMatrix[i][j] = -m_weightedMatrix[j][i];
		}
	}

	while (bellmanFordForFlow(startVertex, weightedMatrix, currentPath))
	{
		int minWeight = INT_MAX;
		currentPath.push_back(m_vertexQuantity - 1);
		int buffer = 0;
		for (int i = 1; i < currentPath.size() / 2; i++)
		{
			buffer = currentPath[i];
			currentPath[i] = currentPath[currentPath.size() - i - 1];
			currentPath[currentPath.size() - i - 1] = buffer;
		}
		//Find maximum flow in path
		for (int i = 0; i < currentPath.size() - 1; i++)
		{
			if (bandwidthMatrix[currentPath[i]][currentPath[i + 1]] < minWeight)
				minWeight = bandwidthMatrix[currentPath[i]][currentPath[i + 1]];
		}
		if (streamValue + minWeight >= streamSize)
		{
			minWeight = streamSize - streamValue;
			isOver = true;
		}
		else streamValue += minWeight;
		//Change matrixes
		for (int i = 0; i < currentPath.size() - 1; i++)
		{
			bandwidthMatrix[currentPath[i]][currentPath[i + 1]] -= minWeight;
			if (bandwidthMatrix[currentPath[i]][currentPath[i + 1]] > 0)
				weightedMatrix[currentPath[i]][currentPath[i + 1]] = modifiedWeightedMatrix[currentPath[i]][currentPath[i + 1]];
			else 
				weightedMatrix[currentPath[i]][currentPath[i + 1]] = 0;

			bandwidthMatrix[currentPath[i + 1]][currentPath[i]] += minWeight;
			if (bandwidthMatrix[currentPath[i + 1]][currentPath[i]] > 0)
				weightedMatrix[currentPath[i + 1]][currentPath[i]] = modifiedWeightedMatrix[currentPath[i + 1]][currentPath[i]];
			else 
				weightedMatrix[currentPath[i + 1]][currentPath[i]] = 0;

			if (currentPath[i] < currentPath[i + 1])
				streamMatrix[currentPath[i]][currentPath[i + 1]] += minWeight;
			if ((streamMatrix[currentPath[i + 1]][currentPath[i]] - minWeight) >= 0)
				streamMatrix[currentPath[i + 1]][currentPath[i]] -= minWeight;
		}
		
		if (isOver)
		{
			std::cout << "Minimal cost flow results: ";
			int count = 0;

			for (int i = 0; i < m_vertexQuantity; i++)
			{
				count += streamMatrix[i][m_vertexQuantity - 1];
			}
			std::cout << "\nFlow: " << count << '\n';
			showMatrix(streamMatrix);
			std::cout << '\n';
			std::cout << "\nCost Matrix: ";
			showMatrix(m_weightedMatrix);
			std::cout << '\n';
			int resultSum = 0;
			for (int i = 0; i < m_vertexQuantity; i++)
			{
				if (streamMatrix[i][m_vertexQuantity - 1] != 0)
				{
					resultSum += streamMatrix[i][m_vertexQuantity - 1] * modifiedWeightedMatrix[i][m_vertexQuantity - 1];
				}
			}
			std::cout << "\nResulting sum: " << resultSum << "\n\n";
			return;
		}
		//Find next path
		currentPath.clear();
		currentPath.push_back(startVertex);
	}

	//Output stream values
	std::cout << "Minimal cost flow results: ";
	int count = 0;

	for (int i = 0; i < m_vertexQuantity; i++)
	{
		count += streamMatrix[i][m_vertexQuantity - 1];
	}
	std::cout << "\nFlow: " << count << '\n';
	showMatrix(streamMatrix);
	std::cout << '\n';
	std::cout << "\nCost Matrix: ";
	showMatrix(m_weightedMatrix);
	std::cout << '\n';
	int resultSum = 0;
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		if (streamMatrix[i][m_vertexQuantity - 1] != 0)
		{
			resultSum += streamMatrix[i][m_vertexQuantity - 1] * modifiedWeightedMatrix[i][m_vertexQuantity - 1];
		}
	}
	std::cout << "\nResulting sum: " << resultSum << "\n\n";
	return;
}

int Graph::prim()
{
	std::mt19937 mersenne(static_cast<unsigned int>(time(0)));
	int iterationCounter = 0;
	m_minimumSpanningTree.clear();
	m_minimumSpanningTree = std::vector<std::vector<int>>(m_vertexQuantity, std::vector<int>(m_vertexQuantity, 0));

	int startingVertex = mersenne() % m_vertexQuantity;
	std::vector<bool> isInSpanningTree(m_vertexQuantity, false);
	isInSpanningTree[startingVertex] = true;
	bool isOver = false;
	
	auto copyOfWeightedMatrix = m_weightedMatrix;
	while(!isOver)
	{
		//Find Minimum Edge connected to current current subgraph
		std::pair<int, int> minimumEdge;
		bool isFound = false;
		for (int i = 0; i < m_vertexQuantity; i++)
		{
			// Check for edges in minimum spanning tree
			if (isInSpanningTree[i] == true)
			{
				for (int j = 0; j < m_vertexQuantity; j++)
				{
					iterationCounter++;
					// Check for not internal edges
					if (isInSpanningTree[j] == false)
					{
						// Check for both possible connections
						if (copyOfWeightedMatrix[i][j] != 0)
						{
							if (!isFound)
							{
								isFound = true;
								minimumEdge = std::make_pair(i, j);
							}
							else
							{
								if (copyOfWeightedMatrix[i][j] < copyOfWeightedMatrix[minimumEdge.first][minimumEdge.second])
								{
									minimumEdge = std::make_pair(i, j);
								}
							}
						}
						if (copyOfWeightedMatrix[j][i] != 0)
						{
							if (!isFound)
							{
								isFound = true;
								minimumEdge = std::make_pair(j, i);
							}
							else
							{
								if (copyOfWeightedMatrix[j][i] < copyOfWeightedMatrix[minimumEdge.first][minimumEdge.second])
								{
									minimumEdge = std::make_pair(j, i);
								}
							}
						}
					}
				}
			}
		}

		//Put edge in minimum spanning tree or not put
		if (!isAchievable(minimumEdge.first, minimumEdge.second, m_minimumSpanningTree, -1))
		{
			m_minimumSpanningTree[minimumEdge.first][minimumEdge.second] = m_weightedMatrix[minimumEdge.first][minimumEdge.second];
			m_minimumSpanningTree[minimumEdge.second][minimumEdge.first] = m_weightedMatrix[minimumEdge.first][minimumEdge.second];
		}
		if (isAchievable(minimumEdge.first, minimumEdge.first, m_minimumSpanningTree, -1))
		{
			m_minimumSpanningTree[minimumEdge.first][minimumEdge.second] = 0;
			m_minimumSpanningTree[minimumEdge.second][minimumEdge.first] = 0;
			copyOfWeightedMatrix[minimumEdge.first][minimumEdge.second] = 0;
		}
		else
		{
			isInSpanningTree[minimumEdge.first] = true;
			isInSpanningTree[minimumEdge.second] = true;
		}

		//Check for end
		isOver = true;
		for (int i = 0; i < m_vertexQuantity; i++)
		{
			if (isInSpanningTree[i] == false)
			{
				isOver = false;
				break;
			}
		}
	}

	std::cout << "\nPrim's Algorithm results:";
	int resultSum = 0;
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		for (int j = i + 1; j < m_vertexQuantity; j++)
		{
			if (m_minimumSpanningTree[i][j] != 0)
			{
				resultSum += m_minimumSpanningTree[i][j];
			}
		}
	}
	std::cout << "\nMinumum Spanning Tree Weight:" << resultSum << '\n';
	showMatrix(m_minimumSpanningTree);
	return iterationCounter;
}

int Graph::kruskal()
{
	int iterationCounter = 0;
	m_minimumSpanningTree.clear();
	m_minimumSpanningTree = std::vector<std::vector<int>>(m_vertexQuantity, std::vector<int>(m_vertexQuantity, 0));
	std::list<std::pair<int, int>> sortedEdges;

	// Fill list with egdes in ascending order
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		for (int j = 0; j < m_vertexQuantity; j++)
		{
			if (m_weightedMatrix[i][j] != 0)
			{
				bool isEmplaced = false;
				for (auto iter = sortedEdges.begin(); iter != sortedEdges.end(); iter++)
				{
					iterationCounter++;
					if (m_weightedMatrix[(*iter).first][(*iter).second] >= m_weightedMatrix[i][j])
					{
						sortedEdges.emplace(iter, std::make_pair(i, j));
						isEmplaced = true;
						break;
					}
				}
				if (!isEmplaced) sortedEdges.push_back(std::make_pair(i, j));
			}
		}
	}

	// Fill minimum spanning tree
	while (!sortedEdges.empty())
	{
		auto currentEdge = *(sortedEdges.begin());
		sortedEdges.pop_front();
		if (!isAchievable(currentEdge.first, currentEdge.second, m_minimumSpanningTree, -1))
		{
			m_minimumSpanningTree[currentEdge.first][currentEdge.second] = m_weightedMatrix[currentEdge.first][currentEdge.second];
			m_minimumSpanningTree[currentEdge.second][currentEdge.first] = m_weightedMatrix[currentEdge.first][currentEdge.second];
		}
		if (isAchievable(currentEdge.first, currentEdge.first, m_minimumSpanningTree, -1))
		{
			m_minimumSpanningTree[currentEdge.first][currentEdge.second] = 0;
			m_minimumSpanningTree[currentEdge.second][currentEdge.first] = 0;
		}
	}

	std::cout << "\nKruskal's Algorithm results:";
	int resultSum = 0;
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		for (int j = i + 1; j < m_vertexQuantity; j++)
		{
			if (m_minimumSpanningTree[i][j] != 0)
			{
				resultSum += m_minimumSpanningTree[i][j];
			}
		}
	}
	std::cout << "\nMinumum Spanning Tree Weight:" << resultSum << '\n';
	showMatrix(m_minimumSpanningTree);
	return iterationCounter;
}

bool Graph::isAchievable(int vertexOne, int vertexTwo, std::vector<std::vector<int>>& graph, int lastVertex)
{
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		if ((graph[i][vertexOne] != 0 || graph[vertexOne][i] != 0) && i != lastVertex)
		{
			if (i == vertexTwo) return true;
			if (isAchievable(i, vertexTwo, graph, vertexOne)) return true;
		}
	}
	return false;
}

int Graph::findNumberOfSpanningTrees()
{
	//Fill Kirchhoff's matrix
	std::vector<std::vector<int>> kirchhoffMatrix(m_vertexQuantity, std::vector<int>(m_vertexQuantity, 0));
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		int degree = 0;
		for (int j = 0; j < m_vertexQuantity; j++)
		{
			if (m_matrix[i][j] != 0)
			{
				kirchhoffMatrix[i][j] = -1;
				degree++;
			}
			if (m_matrix[j][i] != 0)
			{
				kirchhoffMatrix[i][j] = -1;
				degree++;
			}
		}
		kirchhoffMatrix[i][i] = degree;
	}

	//Calculate minor for 1,1
	std::vector<std::vector<int>> newMatrix(kirchhoffMatrix[0].size() - 1, std::vector<int>(kirchhoffMatrix[0].size() - 1, 0));
	for (int j = 0; j < kirchhoffMatrix[0].size() - 1; j++)
	{
		for (int k = 0; k < kirchhoffMatrix[0].size() - 1; k++)
		{
			newMatrix[j][k] = kirchhoffMatrix[j + 1][k + 1];
		}
	}
	return findDeterminant(newMatrix);
}

int Graph::findDeterminant(std::vector<std::vector<int>>& matrix)
{
	if (matrix[0].size() == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	int det = 0;
	for (int i = 0; i < matrix[0].size(); i++)
	{
		std::vector<std::vector<int>> newMatrix(matrix[0].size() - 1, std::vector<int>(matrix[0].size() - 1, 0));
		for (int j = 0; j < matrix[0].size() - 1; j++)
		{
			for (int k = 0; k < matrix[0].size() - 1; k++)
			{
				if (j >= i) newMatrix[j][k] = matrix[j + 1][k + 1];
				else newMatrix[j][k] = matrix[j][k + 1];
			}
		}
		if (i % 2 == 0) det += matrix[0][i] * findDeterminant(newMatrix);
		else det -= matrix[0][i] * findDeterminant(newMatrix);
	}
	return det;
}

void Graph::codePrufer()
{
	m_pruferSpanningTree.clear();

	if (m_minimumSpanningTree.size() == 0)
	{
		std::cout << "\nMinimum spanning tree is not generated.\n";
		return;
	}

	auto copySpanTree = m_minimumSpanningTree;
	int vertexCounter = 0;

	for (int i = 0; i < m_vertexQuantity; i++)
	{
		bool isLeaf = true;
		int count = 0;

		for (int j = 0; j < m_vertexQuantity; j++)
		{
			if (copySpanTree[i][j] != 0)
			{
				count++;
				if (count > 1)
				{
					isLeaf = false;
					break;
				}
			}
		}

		if (isLeaf && count != 0)
		{
			bool isAppripriate = true;
			if (!m_pruferSpanningTree.empty())
			{
				for (auto iter = m_pruferSpanningTree.begin(); iter != m_pruferSpanningTree.begin(); iter++)
				{
					if ((*iter).first == i)
					{
						isAppripriate = false;
						break;
					}
				}
			}

			int nearVertex = 0;
			for (int k = 0; k < m_vertexQuantity; k++)
			{
				if (copySpanTree[i][k] != 0)
				{
					nearVertex = k;
				}
			}

			if (isAppripriate)
			{
				m_pruferSpanningTree.push_back(std::make_pair(i, copySpanTree[i][nearVertex]));

				if (vertexCounter == m_vertexQuantity - 2) 
					m_pruferSpanningTree.push_back(std::make_pair(nearVertex, copySpanTree[i][nearVertex]));

				copySpanTree[i][nearVertex] = 0;
				copySpanTree[nearVertex][i] = 0;
				i = 0;
				vertexCounter++;
			}
		}
	}

	std::cout << "\nPrufer's code: ";
	for (int i = 0; i < m_pruferSpanningTree.size(); i++)
	{
		std::cout << m_pruferSpanningTree[i].first << '\t';
	}
	std::cout << "\n\n";
}

void Graph::decodePrufer()
{
	if (m_pruferSpanningTree.size() == 0)
	{
		std::cout << "\nPrufer's code is not generated.\n";
		return;
	}
}

void Graph::Start()
{
	std::cout << "Matrix without weights:\n";
	showMatrix(m_matrix);
	std::cout << '\n';
	std::cout << "Matrix with weights:\n";
	showMatrix(m_weightedMatrix);
	std::cout << '\n';

	while (true)
	{
		std::cout << "0. Show matrix without weights.\n";
		std::cout << "1. Show table of weighted edges.\n";
		std::cout << "2. Run Shimbell Method.\n";
		std::cout << "3. Show quantity of paths between two vertexes.\n";
		std::cout << "4. Run Dijkstra's Algorithm.\n";
		std::cout << "5. Run Bellman-Ford's Algorithm.\n";
		std::cout << "6. Run Floyd-Warshall's Algorithm.\n";
		std::cout << "7. Show matrix with weights.\n";
		std::cout << "8. Recreate the graph.\n";
		std::cout << "9. Run Ford-Falkerson's Algotithm and Get Flow (2/3 Fmax) of Minimal Cost.\n";
		std::cout << "10. Run Kruskal's Algorithm.\n";
		std::cout << "11. Run Prim's Algorithm.\n";
		std::cout << "12. Find number of Spanning trees using Kirchhoff's theorem.\n";
		std::cout << "13. Code Spanning Tree using Prufer's code.\n";
		std::cout << "14. Decode Prufer's code.\n";
		std::cout << "15. End the program.\n";

		std::string input;
		int operationNumber;
		bool isAssigned = false;

		while (!isAssigned)
		{
			std::cout << "Enter the operation number:\n";
			std::cin >> input;
			if (IsOnlyDigits(input) && std::stoi(input) >= 0 && std::stoi(input) <= 15)
			{
				operationNumber = std::stoi(input);
				isAssigned = true;
			}
			else std::cout << "Number is incorrect.\n";
		}
		switch (operationNumber)
		{
		case 0:
			showMatrix(m_matrix);
			std::cout << '\n';
			break;
		case 1:
			showWeights(m_weightedMatrix);
			std::cout << '\n';
			break;
		case 2:
			{
				int numberOfEdges;
				bool isAssigned = false;
				while (!isAssigned)
				{
					std::cout << "Enter the number of edges in path:\n";
					std::cin >> input;
					if (IsOnlyDigits(input) && std::stoi(input) > 0 && std::stoi(input) < m_vertexQuantity)
					{
						numberOfEdges = std::stoi(input);
						isAssigned = true;
					}
					else std::cout << "Number is incorrect.\n";
				}
				isAssigned = false;
				while (!isAssigned)
				{
					std::cout << "0. Find the longest paths.\n";
					std::cout << "1. Find the shortest paths.\n";
					std::cout << "Enter the number of setting:\n";
					std::cin >> input;
					if (IsOnlyDigits(input) && (std::stoi(input) == 1 || std::stoi(input) == 0))
					{
						shimbellMethod(numberOfEdges, std::stoi(input));
						break;
					}
					else std::cout << "Number is incorrect.\n";
				}
			}
			break;
		case 3:
			{
				int startVertex;
				bool isAssigned = false;
				while (!isAssigned)
				{
					std::cout << "Enter the number of first vertex:\n";
					std::cin >> input;
					if (IsOnlyDigits(input) && std::stoi(input) >= 0 && std::stoi(input) < m_vertexQuantity)
					{
						startVertex = std::stoi(input);
						isAssigned = true;
					}
					else std::cout << "Number is incorrect.\n";
				}
				isAssigned = false;
				while (!isAssigned)
				{
					std::cout << "Enter the number of second vertex:\n";
					std::cin >> input;
					if (IsOnlyDigits(input) && std::stoi(input) >= 0 && std::stoi(input) < m_vertexQuantity)
					{
						std::cout << "Quantity of possible paths: " << dfs(startVertex, std::stoi(input)) << '\n';
						break;
					}
					else std::cout << "Number is incorrect.\n";
				}
			}
			break;
		case 4:
			{
				if (m_mode != 0)
				{
					std::cout << "Dijkstra's algorithm can't be runned in graph with negative weights.\n";
					break;
				}
				bool isAssigned = false;
				while (!isAssigned)
				{
					std::cout << "Enter the number of starting vertex:\n";
					std::cin >> input;
					if (IsOnlyDigits(input) && std::stoi(input) >= 0 && std::stoi(input) < m_vertexQuantity)
					{
						std::cout << "\nNumber of iterations: " << dijkstra(std::stoi(input)) << '\n';
						break;
					}
					else std::cout << "Number is incorrect.\n";
				}
			}
		break;
		case 5:
			{
				bool isAssigned = false;
				while (!isAssigned)
				{
					std::cout << "Enter the number of starting vertex:\n";
					std::cin >> input;
					if (IsOnlyDigits(input) && std::stoi(input) >= 0 && std::stoi(input) < m_vertexQuantity)
					{
						std::cout << "\nNumber of iterations: " << bellmanFord(std::stoi(input)) << '\n';
						break;
					}
					else std::cout << "Number is incorrect.\n";
				}
			}
			break;
		case 6:
			floydWarshall(0);
			break;
		case 7:
			showMatrix(m_weightedMatrix);
			std::cout << '\n';
			break;
		case 8:
			{
				std::string input;
				bool isAssigned = false;
				while (!isAssigned)
				{
					std::cout << "Enter the number of vertexes:\n";
					std::cin >> input;
					if (IsOnlyDigits(input) && std::stoi(input) > 0) isAssigned = true;
					else std::cout << "Number is incorrect.\n";
				}
				createGraph(std::stoi(input));
			}
			std::cout << "Matrix without weights:\n";
			showMatrix(m_matrix);
			std::cout << '\n';
			std::cout << "Matrix with weights:\n";
			showMatrix(m_weightedMatrix);
			std::cout << '\n';
			break;
		case 9:
			{
				bool isAssigned = false;
				while (!isAssigned)
				{
					std::cout << "Enter the number of starting vertex:\n";
					std::cin >> input;
					if (IsOnlyDigits(input) && std::stoi(input) >= 0 && std::stoi(input) < m_vertexQuantity)
					{
						minCostFlow(std::stoi(input), (fordFulkerson(std::stoi(input)) * 2) / 3);
						break;
					}
					else std::cout << "Number is incorrect.\n";
				}
			}
			break;
		case 10:
			{
				int iterations = kruskal();
				std::cout << "Number of iterations: " << iterations << "\n\n";
				break;
			}
		case 11:
			{	
				int iterations = prim();
				std::cout << "Number of iterations: " << iterations << "\n\n";
				break;
			}
		case 12:
			{
				int number = findNumberOfSpanningTrees();
				std::cout << "Number of Spanning Trees: " << number << "\n\n";
				break;
			}
		case 13:
			codePrufer();
			break;
		case 14:
			decodePrufer();
			break;
		case 15:
			return;
		}
	}
}