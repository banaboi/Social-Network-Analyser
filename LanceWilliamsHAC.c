// Implementation of functions defined in LanceWilliamsHAC.h for COMP2521 ass2
// Author : Luke Banicevic, z5209974

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include "LanceWilliamsHAC.h"


#define INF -3.0
#define REMOVED -1
#define CLUSTER -2
#define ALPHA_I .5
#define ALPHA_J .5
#define BETA 0
#define GAMMA .5

////////////////////////////////////////////////////////////////////////
// Function Declerations

static int max(int a, int b);
static double absoluteVal(double a);
static void init_dist(Graph g, double **dist);
static void init_dendA(Graph g, Dendrogram *dendA);
static Dendrogram newNode(Vertex v, Dendrogram clusterI, Dendrogram clusterJ);
static int dendASize(Graph g, Dendrogram *dendA);
static void closestClusters(Graph g, double **dist, int *clusters, Dendrogram *dendA);
static void dendARemove(Graph g, Dendrogram *dendA, int *clusters, Dendrogram *nodes);
static int dendAInsert(Graph g, Dendrogram *dendA, Dendrogram cluster);

static void saveDist(Graph g, double **dist, int *clusters, double *distancesI, double *distancesJ);
static void updateDist(Graph g, double **dist, int index, double *distancesI, double *distancesJ, int method, Dendrogram *dendA);


////////////////////////////////////////////////////////////////////////
// 


/**
 * Generates  a Dendrogram using the Lance-Williams algorithm (discussed
 * in the spec) for the given graph  g  and  the  specified  method  for
 * agglomerative  clustering. The method can be either SINGLE_LINKAGE or
 * COMPLETE_LINKAGE (you only need to implement these two methods).
 * 
 * The function returns a 'Dendrogram' structure.
 */
Dendrogram LanceWilliamsHAC(Graph g, int method) {

    // Create distance array to store distance from vertex i and j
    // and initalise
    double **dist = malloc(sizeof(double *) * GraphNumVertices(g));
    for (int i = 0; i < GraphNumVertices(g); i++) {
        dist[i] = malloc(sizeof(double) * GraphNumVertices(g));
    }
    init_dist(g, dist);

    // Create array where each cell stores a pointer to 1 dendogram which
    // is inialised to store a node for each vertex 
    Dendrogram *dendA = malloc(sizeof(Dendrogram) * GraphNumVertices(g));
    init_dendA(g, dendA);

    
    // While the entire graph hasnt been clustered
    while (dendASize(g, dendA) != 1) {

        // Find closest clusters
        int clusters[2] = {0};
        closestClusters(g, dist, clusters, dendA);
        

        // Remove them from dendA
        Dendrogram nodes[2] = {NULL};
        dendARemove(g, dendA, clusters, nodes);


        // Merge into new cluster
        Dendrogram merged = newNode(CLUSTER, nodes[0], nodes[1]);

        // Add new merged cluster into dendA and return 
        // index where merged cluster will sit in dist array
        int index = dendAInsert(g, dendA, merged);

        // Update distance array by removing closest clusters and save
        // distances for Formula
        double *distancesI = malloc(sizeof(double) * GraphNumVertices(g));
        double *distancesJ = malloc(sizeof(double) * GraphNumVertices(g));
        saveDist(g, dist, clusters, distancesI, distancesJ);
        updateDist(g, dist, index, distancesI, distancesJ, method, dendA);

        // Free memory used this iteration
        free(distancesI);
        free(distancesJ);

    }


    // return final Dendogram
    Dendrogram final = NULL;
    for (int i = 0; i < GraphNumVertices(g); i++) {
        if (dendA[i]->vertex == CLUSTER) {
            final = dendA[i];
            break;
        }
    }
    // Free dist array
    for (int i = 0; i < GraphNumVertices(g); i++) {
        free(dist[i]);
    }
    free(dist);

    // Free dendrogram structure
    for (int i = 0; i < GraphNumVertices(g); i++) {
        if (dendA[i] != final) {
            free(dendA[i]);
        }
    }

    free(dendA);

    return final;
}

// Frees all memory associated with the given Dendrogram structure
void freeDendrogram(Dendrogram d) {
    
    if (d != NULL) {
      freeDendrogram(d->left);
      freeDendrogram(d->right);
      free(d);
   }
}

////////////////////////////////////////////////////////////////////////
// Helper Functions

// Creates new Dendogram Cluster
static Dendrogram newNode(Vertex v, Dendrogram clusterI, Dendrogram clusterJ) {
    Dendrogram newNode = malloc(sizeof(*newNode));

    newNode->vertex = v;
    newNode->left = clusterI;
    newNode->right = clusterJ;

    return newNode;
}

// Initialises the distance array to be used in the Lance Williams Algorithm
static void init_dist(Graph g, double **dist) {
    for (Vertex i = 0; i < GraphNumVertices(g); i++) {
        for (Vertex j = 0; j < GraphNumVertices(g); j++) {
            // If same node, set distance to zero
             if (i == j) {
                dist[i][j] = 0.0;

            // If there is no edge between them, set distance to inf
            }else if (!GraphIsAdjacent(g,i,j) && !GraphIsAdjacent(g,j,i)) {
                dist[i][j] = __DBL_MAX__;
                dist[j][i] = __DBL_MAX__;
            } else {

                int weight = 0;
                // Create list of outs from Vertex i and outs from Vertex j
                AdjList iOut = GraphOutIncident(g, i);
                AdjList jOut = GraphOutIncident(g, j);

                // For all outs, if the dest is j, update max weight
                for (AdjList current = iOut; current != NULL; current = current->next) {
                    if (current->v == j) {
                        weight = max(current->weight, weight);
                    }
                }
                // For all outs, if the dest is i, update max weight
                for (AdjList current = jOut; current != NULL; current = current->next) {
                    if (current->v == i) {
                        weight = max(current->weight, weight);
                    }
                }

                dist[i][j] = (double) 1 / weight;
            }
        }
    }
}

// Initalises the Dendogram to be used in the Lance Williams Algorithm
static void init_dendA(Graph g, Dendrogram *dendA) {

    for (int i = 0; i < GraphNumVertices(g); i++) {
        dendA[i] = newNode(i, NULL, NULL);
    }
}

// Returns the maximum of two numbers
static int max(int a, int b) {
    return (a > b) ? a : b;
}

// Returns the current size of the dendA
static int dendASize(Graph g, Dendrogram *dendA) {

    int size = 0;
    for (int i = 0; i < GraphNumVertices(g); i++) {
        if (dendA[i]->vertex != REMOVED) {
            size++;
        } 
    }

    return size;
}

// Removes the two current closest clusters from the dendA
// and returns them in a nodeClusters array
static void dendARemove(Graph g, Dendrogram *dendA, int *clusters, Dendrogram *nodes) {
        
    
    nodes[0] = dendA[clusters[0]];
    dendA[clusters[0]] = newNode(REMOVED, NULL, NULL);

    nodes[1] = dendA[clusters[1]];
    dendA[clusters[1]] = newNode(REMOVED, NULL, NULL);

}

// Finds the current closest clusters
static void closestClusters(Graph g, double **dist, int *clusters, Dendrogram *dendA) {

    // clusters[0] is i, clusters[1] is j
    double min = __DBL_MAX__;
    for (int i = 0; i < GraphNumVertices(g); i++) {
        for (int j = 0; j < GraphNumVertices(g); j++) {
            // If the current cluster being looked at has been removed, skip
            if (dendA[i]->vertex == REMOVED || dendA[j]->vertex == REMOVED || i == j) continue;
            // if the current distance is smaller than min, update i and j and update
            // currrent min
            if (dist[i][j] < min) {
                min = dist[i][j];
                clusters[0] = i;
                clusters[1] = j;
            }
        }
    }
}

// Adds a new cluster into dendA at the first available location
// and returns this location in the array
static int dendAInsert(Graph g, Dendrogram *dendA, Dendrogram cluster) {

    int index = 0;

    for(int i = 0; i < GraphNumVertices(g); i++) {
        // If free spot, insert
        if (dendA[i]->vertex == REMOVED) {
            dendA[i] = cluster;
            index = i;
            break;
        }
    }

    return index;
}
 
// Saves two arrays of distances to every other cluster from two clusters to calculate new
// distance using 
static void saveDist(Graph g, double **dist, int *clusters, double *distancesI, double *distancesJ) {

    // Remove entries from dist array
    for (int i = 0; i < GraphNumVertices(g); i++) {
        distancesJ[i] = dist[i][clusters[1]];
        
    }
    for (int j = 0; j < GraphNumVertices(g); j++) {
        distancesI[j] = dist[j][clusters[0]];
        
    }
}

// Returns the absolute value of a
static double absoluteVal(double a) {
    return (a < 0) ? a * -1 : a;
}



// Updates the dist array given the newly entered cluster
static void updateDist(Graph g, double **dist, int index, double *distancesI, double *distancesJ, int method, Dendrogram *dendA) {

    for (int k = 0; k < GraphNumVertices(g); k++) {
        // If current cluster being looked at is removed, skip
        if (dendA[k]->vertex == REMOVED) continue;
        // If the distance from I to k or J to k is inf, skip
        if (distancesI[k] != __DBL_MAX__ && distancesJ[k] != __DBL_MAX__) {
            continue;
            
        // If same cluster, distance to self is zero
        } else if (k == index) {
            dist[index][k] = 0.0;

        } else {
            // Calculate Dist(Cij, Ck) using formula provided in spec and update
            if (method == SINGLE_LINKAGE) {
                dist[index][k] = ALPHA_I * distancesI[k] + ALPHA_J * distancesJ[k] - GAMMA * absoluteVal(distancesI[k] - distancesJ[k]);
                dist[k][index] = ALPHA_I * distancesI[k] + ALPHA_J * distancesJ[k] - GAMMA * absoluteVal(distancesI[k] - distancesJ[k]);
            } else {
                dist[index][k] = ALPHA_I * distancesI[k] + ALPHA_J * distancesJ[k] + GAMMA * absoluteVal(distancesI[k] - distancesJ[k]);
                dist[k][index] = ALPHA_I * distancesI[k] + ALPHA_J * distancesJ[k] + GAMMA * absoluteVal(distancesI[k] - distancesJ[k]);
            }
        }
    }
}