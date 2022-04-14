#include "Graph.h"
#include "Constants.h"

#include <iostream>
#include <random>
#include <ctime>
#include <set>

Graph::Graph()
{
	std::cin >> vertexQuantity;
	//Add checks for input

	//Calculate value of p0
	double p0 = 1;
	for (int i = 0; i < vertexQuantity; i++)
	{
		p0 *= r + (i * c);
		p0 /= b + r + (i * c);
	}

	//Calculating Outdegrees
	std::mt19937 mersenne(static_cast<unsigned int>(time(0)));
	for (int i = 0; i < vertexQuantity; i++)
	{
		int x = 0;
		double p = p0;
		double random = static_cast<double>(mersenne() % 10001) / 10000;
		while (true)
		{
			random -= p;
			if (random < 0)
			{
				outdegrees.push_back(x);
				break;
			}
			else
			{
				if (x < vertexQuantity - 1)
				{
					x++;
					p = p * (b + c * (x - 1)) * (vertexQuantity + 1 - x);
					p = p / (random + (vertexQuantity - x) * c) / x;
				}
				else
				{
					outdegrees.push_back(vertexQuantity);
					break;
				}
			}
		}
	}

	//Calculate actual matrix
	for (int i = 0; i < vertexQuantity; i++)
	{
		std::vector<int> temp;
		//Fill with zeros values up to diagonal
		for (int j = 0; j < i + 1; j++) temp.push_back(0);

		//If Outdegree is more than empty space, fill with ones
		if (outdegrees.at(i) >= vertexQuantity - i - 1)
		{
			for (int j = i + 1; j < vertexQuantity; j++) temp.push_back(1);
		}
		else
		{
			std::set<int> pos;
			//Create set with positions of ones
			while (pos.size() != outdegrees.at(i))
			{
				pos.insert((mersenne() % (vertexQuantity - i - 1)) + i + 1);
			}
			//Fill empty space with appropriate values
			for (int j = i + 1; j < vertexQuantity; j++)
			{
				if (pos.find(j) != pos.end()) temp.push_back(1);
				else temp.push_back(0);
			}
		}
		matrix.push_back(temp);
	}
	matrix.at(vertexQuantity - 2).at(vertexQuantity - 1) = 1;
}

void Graph::ShowMatrix()
{
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix.at(i).size(); j++)
		{
			std::cout << matrix.at(i).at(j) << ' ';
		}
		std::cout << '\n';
	}
}