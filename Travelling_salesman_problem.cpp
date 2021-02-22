
#include "stdafx.hpp"

constexpr int V = 4;

// implementation of traveling Salesman Problem
//https://www.geeksforgeeks.org/traveling-salesman-problem-tsp-implementation/
int travellingSalesmanProblem(int graph[][V], int s)
{
	// store all vertex apart from source vertex
	std::vector<int> vertex;
	for (int i = 0; i < V; i++)
		if (i != s)
			vertex.push_back(i);

	// store minimum weight Hamiltonian Cycle.
	int min_path = INT_MAX;
	do {

		// store current Path weight(cost)
		int current_pathweight = 0;

		// compute current path weight
		int k = s;
		for (size_t i = 0; i < vertex.size(); i++) {
			current_pathweight += graph[k][vertex[i]];
			k = vertex[i];
		}
		current_pathweight += graph[k][s];

		// update minimum
		min_path = std::min(min_path, current_pathweight);

	} while (
		std::next_permutation(vertex.begin(), vertex.end()));

	return min_path;
}

// Driver Code
void travellingSalesmanProblem_driver()
{
	// matrix representation of graph
	int graph[][V] = { { 0, 10, 15, 20 },
										 { 10, 0, 35, 25 },
										 { 15, 35, 0, 30 },
										 { 20, 25, 30, 0 } };
	int s = 0;
	std::cout << travellingSalesmanProblem(graph, s) << std::endl;
}