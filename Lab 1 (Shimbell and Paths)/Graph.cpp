#include "Graph.h"
#include "Constants.h"

#include <iostream>
#include <algorithm>
#include <random>
#include <ctime>
#include <set>

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

void Graph::shimbellMethod(int edgeQuantity, bool mode)
{
	auto shimbellMatrix = m_matrix;
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
					if (shimbellMatrix.at(j).at(l) != 0 && shimbellMatrix.at(l).at(k) != 0)
					{
						if (currentValue == 0) currentValue = shimbellMatrix.at(j).at(l) + shimbellMatrix.at(l).at(k);
						else
						{
							if (mode) currentValue = std::min(currentValue, shimbellMatrix.at(j).at(l) + shimbellMatrix.at(l).at(k));
							else currentValue = std::max(currentValue, shimbellMatrix.at(j).at(l) + shimbellMatrix.at(l).at(k));
						}
					}
				}
				tempMatrix.at(j).at(k) = currentValue;
			}
		}
		shimbellMatrix = tempMatrix;
	}

	showMatrix(shimbellMatrix);
}

void Graph::Start()
{
	// Add Checks for shimbell input
	showMatrix(m_matrix);
	shimbellMethod(2, true);
	shimbellMethod(4, true);
}