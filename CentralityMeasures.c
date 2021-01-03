// Implementation of functions defined in CentralityMeasures.h for COMP2521 ass2
// Author : Luke Banicevic, z5209974

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "CentralityMeasures.h"
#include "Dijkstra.h"

////////////////////////////////////////////////////////////////////////
// Function Declerations

static NodeValues newNodeValues(Graph g);
static double shortestPathDistance(ShortestPaths paths, Vertex u);
static double nodesReachable(ShortestPaths paths);
static double calculateClosenessCentrality(ShortestPaths paths);
static double calculateBetCentrality(Graph g, ShortestPaths paths);
static double calculateNormBetCentrality(Graph g, double betCentrality);
static void countPaths(ShortestPaths paths, PredNode *node, Vertex v, double *noPaths, double *pathsV);

////////////////////////////////////////////////////////////////////////



// Finds the closeness centrality for each vertex in the given graph and
// returns the results in a NodeValues structure.
NodeValues closenessCentrality(Graph g) {

    // Create the NodeValues structure
    NodeValues n = newNodeValues(g);
    
    // Loop through all vertices in graph
    for (Vertex u = 0; u < GraphNumVertices(g); u++) {

        // For each vertice, create a shortestPath structure 
        ShortestPaths paths = dijkstra(g, u);

        // Fill closeness centrality value using paths structure
        n.values[u] = calculateClosenessCentrality(paths);
    }

    return n;

}


// Finds  the  betweenness centrality for each vertex in the given graph
// and returns the results in a NodeValues structure.
NodeValues betweennessCentrality(Graph g) {

    // Create the NodeValues structure
    NodeValues n = newNodeValues(g);

    for (Vertex u = 0; u < GraphNumVertices(g); u++) {

        // For each vertice, create a shortestPath structure 
        ShortestPaths paths = dijkstra(g, u);
        
        // Fill betweenness centrality value using paths structure
        n.values[u] = calculateBetCentrality(g, paths);

    }

    return n;
}


// Finds  the  normalised  betweenness centrality for each vertex in the
// given graph and returns the results in a NodeValues structure.
NodeValues betweennessCentralityNormalised(Graph g) {
    // Create the NodeValues structure
    NodeValues n = newNodeValues(g);

    for (Vertex u = 0; u < GraphNumVertices(g); u++) {

        // For each vertice, create a shortestPath structure 
        ShortestPaths paths = dijkstra(g, u);
        
        // Fill norm betweenness centrality value using paths structure
        double betCentrality = calculateBetCentrality(g, paths);
        n.values[u] = calculateNormBetCentrality(g, betCentrality);

    }

    return n;
}

// Prints out the nodevalues for debugging
void showNodeValues(NodeValues nvs) {
    for (int i = 0; i < nvs.numNodes; i++) {
        printf("%d: %f\n",i, nvs.values[i]);
    }
}

// Frees all memory associated with the given NodeValues structure.
void freeNodeValues(NodeValues nvs) {
    free(nvs.values);
}

////////////////////////////////////////////////////////////////////////
// Helper functions


// Creates a NodeValues structure
static NodeValues newNodeValues(Graph g) {

    NodeValues newNodeValues;

    newNodeValues.numNodes = GraphNumVertices(g);
    newNodeValues.values = malloc(sizeof(double) * GraphNumVertices(g));
    for (Vertex i = 0; i < newNodeValues.numNodes; i++) {
        newNodeValues.values[i] = 0.0;
    }

    return newNodeValues;
}

// Calculates the shortest-path distance in a directed graph from 
// vertex u to v
static double shortestPathDistance(ShortestPaths paths, Vertex u) {
    
    double distance = 0.0;
    for (Vertex v = 0; v < paths.numNodes; v++) {
        // If v is reachable from u, add distance
        if (paths.dist[v] != 0) {
            distance += paths.dist[v];
        }
    }


    return distance;
}

// Calculates the closeness centrality of a node u (paths.src) in a directed graph
static double calculateClosenessCentrality(ShortestPaths paths) {

    // Nodes reachable from u
    double n = nodesReachable(paths);

    // If there are no nodes reachable, return 0.0;
    if (!n) return 0.0;

    // Number of nodes in graph
    double N = paths.numNodes;

    // Shortest-path distance sum of all v reachable from u
    double distance = shortestPathDistance(paths, paths.src);

    // If the sum of all v reachable from u is zero, return zero
    if (distance == 0.0) return 0.0;


    return ((n - 1) / (N - 1)) * ((n - 1) / distance);
    
}

// Calculates the betweenness centrality of a node v (paths.src) in a directed graph
/*  Requirements:
    - total number of shortest paths from s to t for every pair of s and t in g->nV, where s != t != v
    - number of paths which pass through v
*/

static double calculateBetCentrality(Graph g, ShortestPaths paths) {
    
    // Define variables
    double noPaths = 0;
    double pathsV = 0;
    double betCentrality = 0;
    Vertex v = paths.src;


    
    // Generate all pairs of s and t in g such that s != t != v
    for (Vertex s = 0; s < paths.numNodes; s++) {
        // If s is v, skip
        if (s == v) continue;
        // Generate shortest path struct for src s
        ShortestPaths paths_s = dijkstra(g, s);
        for (Vertex t = 0; t < paths.numNodes; t++) {
            // If t is v or t is s, skip
            if (t == v || t == s) continue;
            // If there exists a path between s and t, update variables
            // and calculate betweenness centrality
            if (paths_s.dist[t]) {
                countPaths(paths_s, paths_s.pred[t], v, &noPaths, &pathsV);
                // Accumulate sum
                betCentrality += (pathsV / noPaths);
                // Reset variables for next pair of vertices s and t
                pathsV = 0;
                noPaths = 0;
            }
            
            
        }
    }


   return betCentrality;

}
// Calculates the normalised betweeness centrality of a node in a graph
static double calculateNormBetCentrality(Graph g, double betCentrality) {
    return betCentrality / ((GraphNumVertices(g)-1) * (GraphNumVertices(g) - 2));
}


// Returns the number of nodes reachable from u, paths.src
static double nodesReachable(ShortestPaths paths) {

    double n = 0.0;
    for (Vertex v = 0; v < paths.numNodes; v++) {

        // If v is reachable from u, add
        if (paths.dist[v] != 0 || v == paths.src) {
            n++;
        }
    }

    return n;
}

// Recursive function to count the number of shortest paths from s to t,
// including those paths which contain the vertex v
static void countPaths(ShortestPaths paths, PredNode *node, Vertex v, double *noPaths, double *pathsV) {
    
    // Base case: prednode is NULL
    if (!node) {
        return;
        // Case 1: src is reached = increment noPaths
    } else if (node->v == paths.src) {
        *noPaths += 1.0;
        // Case 2: v is visited where v != src, increment pathsV
    } else if (node->v == v && v != paths.src) {
        *pathsV += 1.0;
    }

    // case 3: move through predecessors
    countPaths(paths, paths.pred[node->v], v, noPaths, pathsV);
    countPaths(paths, node->next, v, noPaths, pathsV);
}