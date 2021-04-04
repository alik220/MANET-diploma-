#pragma once

#include "Header.h"
#include <stack>
#include <string>
#include <iostream>

using namespace std;

double betw_cent(int x1, int y1, int x2, int y2, int *playerPos, int quantAgent);
void find_neighbor(int x, int y, int* playerPos, int quantAgent);
//void BrandesAlgorithm_Unweighted(double CB[], NETWORK* network, int quantAgent);

void find_neighbor(int* array, int *playerPos, int quantAgent)
{
	int** net = new int* [quantAgent]; //Воссоздание сети
	for (int i = 0; i < quantAgent; i++)
	{
		net[i] = new int[2];
	}

	for (int i = 0; i < quantAgent; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			net[i][j] = *playerPos;
			playerPos++;
		}
	}

	int vertex[2];
	vertex[0] = *array;
	array++;
	vertex[1] = *array;

	for (int i; i < quantAgent; i++)
	{

	}
}

/*
void BrandesAlgorithm_Unweighted(double CB[], NETWORK* network, int quantAgent) {

	double i, j, u, v;

	vector<unsigned long> d;								// A vector storing shortest distance estimates
	vector<unsigned long> sigma;							// sigma is the number of shortest paths
	vector<double> delta;							// A vector storing dependency of the source vertex on all other vertices
	vector< vector <unsigned long> > PredList;			// A list of predecessors of all vertices 

	queue <unsigned long> Q;								// A priority queue soring vertices
	stack <unsigned long> S;								// A stack containing vertices in the order found by Dijkstra's Algorithm


	// Compute Betweenness Centrality for every vertex i
	for (i = 0; i < quantAgent; i++) {

		
		PredList.assign(quantAgent, vector <unsigned long>(0, 0));
		d.assign(quantAgent, ULONG_MAX);
		d[i] = 0;
		sigma.assign(quantAgent, 0);
		sigma[i] = 1;
		delta.assign(quantAgent, 0);
		Q.push(i);

		// Use Breadth First Search algorithm 
		while (!Q.empty()) {
			// Get the next element in the queue
			u = Q.front();
			Q.pop();
			// Push u onto the stack S. Needed later for betweenness computation
			S.push(u);
			// Iterate over all the neighbors of u 
			for (j = 0; j < (unsigned long)network->vertex[u].degree; j++) {
				// Get the neighbor v of vertex u
				v = (unsigned long)network->vertex[u].edge[j].target;

				
				if (d[v] == ULONG_MAX) {
					d[v] = d[u] + 1;
					Q.push(v);
				}
				if (d[v] == d[u] + 1) {
					sigma[v] += sigma[u];
					PredList[v].push_back(u);
				}
			} // End For

		} // End While 

		
		while (!S.empty()) {
			u = S.top();
			S.pop();
			for (j = 0; j < PredList[u].size(); j++) {
				delta[PredList[u][j]] += ((f64)sigma[PredList[u][j]] / sigma[u]) * (1 + delta[u]);
			}
			if (u != i)
				CB[u] += delta[u];
		}

		// Clear data for the next run
		PredList.clear();
		d.clear();
		sigma.clear();
		delta.clear();

	} // End For 

	// End time after Brandes' algorithm and the time difference
	time(&end);
	time_dif = difftime(end, start);
	cout << "It took " << time_dif << " seconds to calculate Betweenness Centrality in an unweighted graph" << endl;

	return;

} // End of BrandesAlgorithm_Unweighted 

*/

double betw_cent(int x1, int y1, int x2, int y2, int* playerPos, int quantAgent)
{
	int** net2Drone = new int* [quantAgent + 2]; //Воссоздание сети с учетом пары дронов
	for (int i = 0; i < quantAgent + 2; i++)
	{
		net2Drone[i] = new int[2];
	}

	for (int i = 0; i < quantAgent; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			net2Drone[i][j] = *playerPos;
			playerPos++;
		}
	}

	net2Drone[quantAgent][0] = x1;
	net2Drone[quantAgent + 1][0] = x2;
	net2Drone[quantAgent][1] = y1;
	net2Drone[quantAgent + 1][1] = y2;



	return 0;
}
