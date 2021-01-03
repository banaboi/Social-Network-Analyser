
// Implementation of Dijkstra's Algorithm for COMP2521 ass2
// Author : Luke Banicevic, z5209974

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "Dijkstra.h"
#include "PQ.h"

////////////////////////////////////////////////////////////////////////
// Function Declerations

static ShortestPaths newShortestPaths(Graph g, Vertex src);

static PredNode *newPredNode(Vertex v);
static PredNode *predNodeInsert(ShortestPaths paths, PredNode *p, Vertex v);
static void predNodeFree(PredNode *p);
static PredNode *predNodeFilter(Graph g, ShortestPaths paths, PredNode *p, Vertex dest);

////////////////////////////////////////////////////////////////////////

/*
 * The  function  returns  a 'ShortestPaths' structure with the required
 * information:
 * - the number of vertices in the graph
 * - the source vertex
 * - distance array
 * - array of predecessor lists
 */
ShortestPaths dijkstra(Graph g, Vertex src) {

	// Create shortestPath structure
    ShortestPaths paths = newShortestPaths(g, src);

    // If the graph doesnt exist, return NULL
    if (g == NULL) return paths;

	// Distance from src to src is zero
	paths.dist[src] = 0;

	// Create a priority queue 
	PQ pq = PQNew();

	// Place src into pq
	PQInsert(pq, src, paths.dist[src]);

	// While pq isnt empty
	while (!PQIsEmpty(pq)) {
		
		// Extract vertex from pq
		Vertex u = PQDequeue(pq);

		// Create list of adjacent nodes coming out of u
		AdjList uOut = GraphOutIncident(g, u);

		// For all adjacent nodes, update distance and predecessors
		for (AdjList current = uOut; current != NULL; current = current->next) {
			Vertex v = current->v;
			// If the current distance set for v is larger, relax
			if (paths.dist[v] >= paths.dist[u] + current->weight) {
				paths.dist[v] = paths.dist[u] + current->weight;
				PQInsert(pq, v, paths.dist[v]);
				paths.pred[v] = predNodeInsert(paths, paths.pred[v], u);
				//DEBUGGING
				//printf("%d added to pred list of %d\n", u, v);
			}

		}

	}

	// Filter pred list for each node
	for (Vertex v = 0; v < GraphNumVertices(g); v++) {
		paths.pred[v] = predNodeFilter(g, paths, paths.pred[v], v);
	}

	
	// Set all distances to unreachable vertices to zero
	for (Vertex v = 0; v < paths.numNodes; v++) {
		if (paths.dist[v] == __INT_MAX__ && v != src) {
			paths.dist[v] = 0;
		}
	}

	
	PQFree(pq);

	return paths;

}

// Frees all memory associated with the given ShortestPaths structure.
void freeShortestPaths(ShortestPaths sps) {
	// Free distance array
	free(sps.dist);

	// Free pred array of PredNodes
	for(int i = 0; i < sps.numNodes; i++) {
		predNodeFree(sps.pred[i]);
	}

	free(sps.pred);
}

////////////////////////////////////////////////////////////////////////
// Helper functions


// Creates a new ShortestPath 
static ShortestPaths newShortestPaths(Graph g, Vertex src) {
    ShortestPaths newPaths;
    newPaths.numNodes = GraphNumVertices(g);
    newPaths.src = src;
    newPaths.dist = malloc(sizeof(Vertex) * GraphNumVertices(g)); 
    newPaths.pred = malloc(sizeof(PredNode) * GraphNumVertices(g));


	// Intialise pred and distance arrays
	for (int i = 0; i < GraphNumVertices(g); i++) {
		newPaths.dist[i] = __INT_MAX__;
		newPaths.pred[i] = NULL;
	}

    return newPaths;

}

// Creates a new PredNode for a PredNode list
static PredNode *newPredNode(Vertex v) {
	
	PredNode *newPredNode = malloc(sizeof(PredNode));
	newPredNode->v = v;
	newPredNode->next = malloc(sizeof(PredNode));
	newPredNode->next = NULL;

	return newPredNode;
}

// Appends a new PredNode into the given PredNode List
static PredNode *predNodeInsert(ShortestPaths paths, PredNode *p, Vertex v) {
	
	// If list is empty, insert at head
	if (p == NULL) {
		return newPredNode(v);
	}
	// Otherwise append to the list
	PredNode *tail = p;
	while (tail->next) {
		tail = tail->next;
	}

	tail->next = newPredNode(v);
	return p;
}

	
// Frees all memory associated with the given PredNode list
static void predNodeFree(PredNode *p) {
	if (p == NULL) {
		return;
	}

	PredNode *current = p;
	while (current) {
		PredNode *next = current->next;
		free(current);
		current = next;
	}
}
// Removes invalid shortest path predecessors after completion of Dijkstras Algorithm
static PredNode *predNodeFilter(Graph g, ShortestPaths paths, PredNode *p, Vertex dest) {

	if (p == NULL) {
		return NULL;
	}

	int size = 0;
	for (PredNode *current = p; current; current = current->next) {
		size++;
	}

	if (size == 1) return p;

	int *weights = malloc(sizeof(int) * GraphNumVertices(g));
	for (Vertex v = 0; v < GraphNumVertices(g); v++) {
		weights[v] = 0;
	}
	
	// Find min weight into dest
	int min = __INT_MAX__;
	AdjList destIn = GraphInIncident(g, dest);
	for (AdjList current = destIn; current != NULL; current = current->next) {
		weights[current->v] = current->weight;
		if (current->weight + paths.dist[current->v] < min) {
			min = current->weight + paths.dist[current->v];
		}
	}

	

	// Loop through pred list, if the one of the preds has a bigger weight
	// to dest then min, remove
	PredNode *current = p;
	while (current && weights[current->v] + paths.dist[current->v] != min) {
		PredNode *delete = current;
		current = current->next;
		free(delete);
	} 
	
	
	free(weights);
	return current;
}