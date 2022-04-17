#include "Graph.h"
#include "Constants.h"

#include <iostream>
#include <algorithm>
#include <random>
#include <ctime>
#include <set>
#include <queue>
#include <deque>

#define stop __asm nop

Graph::Graph()
{
	std::cin >> m_vertexQuantity;
	//Add checks for input

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
	m_weightedMatrix = m_matrix;
	for (int i = 0; i < m_weightedMatrix.size(); i++)
	{
		for (int j = 0; j < m_weightedMatrix.at(i).size(); j++)
		{
			if (m_weightedMatrix.at(i).at(j) == 1) m_weightedMatrix.at(i).at(j) = (mersenne() % 100) + 1;
		}
	}
}

void Graph::showMatrix(std::vector<std::vector<int>>& matrix)
{
	std::cout << '\n';
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix.at(i).size(); j++)
		{
			std::cout << matrix.at(i).at(j) << ' ';
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

	//Output marks
	std::cout << "Dijkstra results:\n";
	for (int i = 0; i < m_vertexQuantity; i++)
	{
		std::cout << startVertex << "-->" << i << " Shortest path length: ";
		if (marks[i] == INT_MAX) std::cout << "infinity";
		else std::cout << marks[i];
		std::cout << '\n';
	}
	return iterationCounter;
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
		std::cout << startVertex << "-->" << i << " Shortest path length: ";
		if (marks[i] == INT_MAX) std::cout << "infinity";
		else std::cout << marks[i];
		std::cout << '\n';
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
			std::cout << j << "-->" << i << " Shortest path length: ";
			if (distance[i][j] == INT_MAX) std::cout << "infinity";
			else std::cout << distance[i][j];
			std::cout << '\n';
		}
	}
	return iterationCounter;
}

void Graph::Start()
{
	// Add Checks for shimbell input
	showMatrix(m_matrix);
	showWeights(m_weightedMatrix);
	//shimbellMethod(2, true);
	//shimbellMethod(3, true);
	while (true)
	{
	/*	int first;
		int second;
		std::cin >> first >> second;
		if (first == -1) break;
		std::cout << dfs(first, second);*/
		int vertex = 0;
		std::cin >> vertex;
		if (vertex == -1) break;
		std::cout << '\n' << dijkstra(vertex) << '\n';
		std::cout << '\n' << bellmanFord(vertex);
		std::cout << '\n' << floydWarshall(vertex);
	}
}