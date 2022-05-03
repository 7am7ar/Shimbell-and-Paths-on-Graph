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

#define stop __asm nop

Graph::Graph()
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

void Graph::createShit()
{
	m_matrix.clear();
	m_weightedMatrix.clear();
	m_outdegrees.clear();
	m_vertexQuantity = 4;
	std::vector<int> one = { 0, 1000, 1000, 0 };
	std::vector<int> two = { 0, 0, 1, 1000 };
	std::vector<int> three = { 0, 0, 0, 1000 };
	std::vector<int> four = { 0, 0, 0, 0};
	m_weightedMatrix.push_back(one);
	m_weightedMatrix.push_back(two);
	m_weightedMatrix.push_back(three);
	m_weightedMatrix.push_back(four);
}

void Graph::createGraph(int vertexQuantity)
{
	// Delete the previous graph
	m_matrix.clear();
	m_weightedMatrix.clear();
	m_outdegrees.clear();
	m_vertexQuantity = vertexQuantity;

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
}


void Graph::showMatrix(std::vector<std::vector<int>>& matrix)
{
	std::cout << '\n';
	for (int i = 0; i < matrix.size(); i++)
	{
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
		if (marks[i] + m_weightedMatrix[i][finalVertex] == marks[finalVertex])
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

bool Graph::dfsFordFulkerson(std::vector<std::vector<int>>& graph, std::vector<int>& currentPath, int startVertex)
{
	bool result = false;

	for (int i = 0; i < m_vertexQuantity; i++)
	{
		if (graph[startVertex][i] > 0)
		{
			currentPath.push_back(i);
			if (i == m_vertexQuantity - 1) return true;
			result = dfsFordFulkerson(graph, currentPath, i);
			if (result == false) currentPath.pop_back();
			else return true;
		}
	}
	return false;
}

void Graph::fordFulkerson(int startVertex)
{
	std::vector<std::vector<int>> streamMatrix(m_vertexQuantity, std::vector<int>(m_vertexQuantity, 0));
	auto bandwidthMatrix = m_weightedMatrix;
	std::vector<int> currentPath;
	currentPath.push_back(startVertex);

	//Find path or not find
	bool result = dfsFordFulkerson(bandwidthMatrix, currentPath, startVertex);;

	while (result)
	{
		int minWeight = INT_MAX;
		//Find minimum weight in path
		for (int i = 0; i < currentPath.size() - 1; i++)
		{
			if (bandwidthMatrix[currentPath[i]][currentPath[i + 1]] < minWeight)
				minWeight = bandwidthMatrix[currentPath[i]][currentPath[i + 1]];
		}

		//Change matrixes
		for (int i = 0; i < currentPath.size() - 1; i++)
		{
			bandwidthMatrix[currentPath[i]][currentPath[i + 1]] -= minWeight;
			bandwidthMatrix[currentPath[i + 1]][currentPath[i]] -= minWeight;
			streamMatrix[currentPath[i]][currentPath[i + 1]] += minWeight;
		}
		//Find next path
		currentPath.clear();
		currentPath.push_back(startVertex);
		result = dfsFordFulkerson(bandwidthMatrix, currentPath, startVertex);
	}

	//Output stream values
	std::cout << "Ford-Fulkerson results: ";
	showMatrix(streamMatrix);
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
		std::cout << "9. Run Ford-Falkerson's Algotithm.\n";
		std::cout << "10. End the program.\n";

		std::string input;
		int operationNumber;
		bool isAssigned = false;

		while (!isAssigned)
		{
			std::cout << "Enter the operation number:\n";
			std::cin >> input;
			if (IsOnlyDigits(input) && std::stoi(input) >= 0 && std::stoi(input) <= 11)
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
						fordFulkerson(std::stoi(input));
						break;
					}
					else std::cout << "Number is incorrect.\n";
				}
			}
			break;
		case 10:
			return;
		case 11:
			createShit();
			break;
		}
	}
}