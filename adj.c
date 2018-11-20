#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEFAULT_ARR_SIZE 20
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// A structure to represent a node in adjacency list
typedef struct AdjListNode {
  int dest;
  int weight;
  struct AdjListNode* next;
} AdjListNode;

// A structure to represent an adjacency liat
typedef struct AdjList {
  AdjListNode* head;  // pointer to head node of list
} AdjList;

// A structure to represent a graph. A graph is an array of adjacency lists.
// Size of array will be V (number of vertices in graph)
typedef struct Graph {
  int V, E;
  AdjList* array;
} Graph;

// A utility function to create a new adjacency list node
AdjListNode* newAdjListNode(int dest, int weight) {
  AdjListNode* newNode = (AdjListNode*)malloc(sizeof(AdjListNode));
  newNode->dest = dest;
  newNode->weight = weight;
  newNode->next = NULL;
  return newNode;
}

void destroyAdjList(AdjListNode* root) {
  AdjListNode *next = NULL, *current = root;
  while (current != NULL) {
    next = current->next;
    free(current);
    current = next;
  }
}

typedef struct DistNode {
  int dist;
  int from, steps;
} DistNode;

// Adds an edge to an undirected graph
void addEdge(Graph* graph, int src, int dest, int weight) {
  // Add an edge from src to dest.  A new node is added to the adjacency
  // list of src.  The node is added at the begining
  AdjListNode* newNode = newAdjListNode(dest, weight);
  newNode->next = graph->array[src].head;
  graph->array[src].head = newNode;

  // Since graph is undirected, add an edge from dest to src also
  newNode = newAdjListNode(src, weight);
  newNode->next = graph->array[dest].head;
  graph->array[dest].head = newNode;
}

// Structure to represent a min heap node
typedef struct MinHeapNode {
  int v;
  int dist;
} MinHeapNode;

// Structure to represent a min heap
typedef struct MinHeap {
  int size;      // Number of heap nodes present currently
  int capacity;  // Capacity of min heap
  int* pos;      // This is needed for decreaseKey()
  MinHeapNode** array;
} MinHeap;

// A utility function to create a new Min Heap Node
MinHeapNode* newMinHeapNode(int v, int dist) {
  MinHeapNode* minHeapNode = (MinHeapNode*)malloc(sizeof(MinHeapNode));
  minHeapNode->v = v;
  minHeapNode->dist = dist;
  return minHeapNode;
}

// A utility function to create a Min Heap
MinHeap* createMinHeap(int capacity) {
  MinHeap* minHeap = (MinHeap*)malloc(sizeof(MinHeap));
  minHeap->pos = (int*)malloc(capacity * sizeof(int));
  minHeap->size = 0;
  minHeap->capacity = capacity;
  minHeap->array = (MinHeapNode**)malloc(capacity * sizeof(MinHeapNode*));
  return minHeap;
}

void cleanMinHeap(MinHeap* heap) {
  int i = 0;
  for (i = 0; i < heap->size; i++) {
    free(heap->array[i]);
  }

  for (i = 0; i < heap->capacity; i++) {
    heap->pos[i] = 0;
  }
  heap->size = 0;
}

void destroyMinHeap(MinHeap* heap) {
  int i = 0;
  for (i = 0; i < heap->size; i++) {
    free(heap->array[i]);
  }
  free(heap->array);
  free(heap->pos);
  free(heap);
}

// A utility function to swap two nodes of min heap. Needed for min heapify
void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b) {
  MinHeapNode* t = *a;
  *a = *b;
  *b = t;
}

// A standard function to heapify at given idx
// This function also updates position of nodes when they are swapped.
// Position is needed for decreaseKey()
void minHeapify(MinHeap* minHeap, int idx) {
  int smallest, left, right;
  smallest = idx;
  left = 2 * idx + 1;
  right = 2 * idx + 2;

  if (left < minHeap->size &&
      minHeap->array[left]->dist < minHeap->array[smallest]->dist)
    smallest = left;

  if (right < minHeap->size &&
      minHeap->array[right]->dist < minHeap->array[smallest]->dist)
    smallest = right;

  if (smallest != idx) {
    // The nodes to be swapped in min heap
    MinHeapNode* smallestNode = minHeap->array[smallest];
    MinHeapNode* idxNode = minHeap->array[idx];

    // Swap positions
    minHeap->pos[smallestNode->v] = idx;
    minHeap->pos[idxNode->v] = smallest;

    // Swap nodes
    swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);

    minHeapify(minHeap, smallest);
  }
}

// A utility function to check if the given minHeap is ampty or not
int isEmpty(MinHeap* minHeap) { return minHeap->size == 0; }

// Standard function to extract minimum node from heap
MinHeapNode* extractMin(MinHeap* minHeap) {
  if (isEmpty(minHeap)) return NULL;

  // Store the root node
  MinHeapNode* root = minHeap->array[0];

  // Replace root node with last node
  MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
  minHeap->array[0] = lastNode;

  // Update position of last node
  minHeap->pos[root->v] = minHeap->size - 1;
  minHeap->pos[lastNode->v] = 0;

  // Reduce heap size and heapify root
  --minHeap->size;
  minHeapify(minHeap, 0);

  return root;
}

// Function to decreasy dist value of a given vertex v. This function
// uses pos[] of min heap to get the current index of node in min heap
void decreaseKey(MinHeap* minHeap, int v, int dist) {
  // Get the index of v in  heap array
  int i = minHeap->pos[v];

  // Get the node and update its dist value
  minHeap->array[i]->dist = dist;

  // Travel up while the complete tree is not hepified.
  // This is a O(Logn) loop
  while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist) {
    // Swap this node with its parent
    minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
    minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
    swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);

    // move to parent index
    i = (i - 1) / 2;
  }
}

// A utility function to check if a given vertex
// 'v' is in min heap or not
bool isInMinHeap(MinHeap* minHeap, int v) {
  if (minHeap->pos[v] < minHeap->size) return true;
  return false;
}

// The main function that calulates distances of shortest paths from src to all
// vertices. It is a O(ELogV) function
int dijkstra(struct Graph* graph, int src, int dest, DistNode* distArr,
             MinHeap* minHeap) {
  int V = graph->V;  // Get the number of vertices in graph

  // minHeap represents set E

  // Initialize min heap with all vertices. dist value of all vertices
  int v = 0;
  for (v = 0; v < V; ++v) {
    distArr[v].dist = INT_MAX;
    distArr[v].from = src;
    distArr[v].steps = INT_MAX;
    minHeap->array[v] = newMinHeapNode(v, distArr[v].dist);
    minHeap->pos[v] = v;
  }

  // Make dist value of src vertex as 0 so that it is extracted first
  minHeap->array[src] = newMinHeapNode(src, distArr[src].dist);
  minHeap->pos[src] = src;
  distArr[src].dist = 0;
  distArr[src].steps = 1;
  distArr[src].from = src;
  decreaseKey(minHeap, src, distArr[src].dist);

  // Initially size of min heap is equal to V
  minHeap->size = V;

  // In the followin loop, min heap contains all nodes
  // whose shortest distance is not yet finalized.
  while (!isEmpty(minHeap)) {
    // Extract the vertex with minimum distance value
    MinHeapNode* minHeapNode = extractMin(minHeap);
    int u = minHeapNode->v;  // Store the extracted vertex number
    free(minHeapNode);

    // Traverse through all adjacent vertices of u (the extracted
    // vertex) and update their distance values
    AdjListNode* pCrawl = graph->array[u].head;
    while (pCrawl != NULL) {
      int v = pCrawl->dest;

      // If shortest distance to v is not finalized yet, and distance to v
      // through u is less than its previously calculated distance
      if (isInMinHeap(minHeap, v) && distArr[u].dist != INT_MAX &&
          pCrawl->weight + distArr[u].dist < distArr[v].dist) {
        distArr[v].dist = distArr[u].dist + pCrawl->weight;
        distArr[v].from = u;
        distArr[v].steps = distArr[u].steps + 1;

        // update distance value in min heap also
        decreaseKey(minHeap, v, distArr[v].dist);
      }
      pCrawl = pCrawl->next;
    }
  }

  // print the calculated shortest distances
  return distArr[dest].dist;
}

Graph* createGraph(int V, int E) {
  Graph* graph = (Graph*)malloc(sizeof(Graph));

  graph->V = V;
  graph->E = E;
  graph->array = (AdjList*)malloc(V * sizeof(AdjList));

  int i;
  for (i = 0; i < V; ++i) {
    graph->array[i].head = NULL;
  }

  return graph;
}

void destroyGraph(Graph* graph) {
  int i;
  AdjListNode* temp = NULL;
  for (i = 0; i < graph->V; ++i) {
    while (graph->array[i].head != NULL) {
      temp = graph->array[i].head;
      graph->array[i].head = temp->next;
      free(temp);
    }
  }

  // destroyAdjList(graph->array->head);
  free(graph->array);
  free(graph);
}

typedef struct Node {
  int x, y;
} Node;

Graph* read_graph(char* filename, int* numVertex, int* numEdge) {
  FILE* fptr;

  fptr = fopen(filename, "r");
  if (fptr == NULL) {
    printf("Cannot open the file %s\n", filename);
    exit(1);
  }

  fscanf(fptr, "%d %d", numVertex, numEdge);
  Graph* g = createGraph(*numVertex, *numEdge);

  Node* nodes = (Node*)malloc((*numVertex) * sizeof(Node));
  int i = 0;
  int j, x, y;
  for (i = 0; i < *numVertex; i++) {
    fscanf(fptr, "%d %d %d", &j, &x, &y);
    nodes[j].x = x;
    nodes[j].y = y;
  }

  int src, dest;
  int _x, _y, weight;

  for (i = 0, _x = 0.0, _y = 0.0; i < *numEdge; i++) {
    fscanf(fptr, "%u %u", &src, &dest);
    _x = abs(nodes[src].x - nodes[dest].x);
    _y = abs(nodes[src].y - nodes[dest].y);

    weight = pow(pow(_x, 2.0) + pow(_y, 2.0), 0.5);

    addEdge(g, src, dest, weight);
  }

  free(nodes);
  fclose(fptr);
  return g;
}

void read_query(char* filename, Graph* graph, DistNode* distArr) {
  FILE* fptr;

  fptr = fopen(filename, "r");
  if (fptr == NULL) {
    printf("Cannot open the file %s\n", filename);
    exit(1);
  }

  // Step 1: Initialize distances from src to all other vertices
  // as INFINITE
  int cases = 0;

  fscanf(fptr, "%d", &cases);
  int i = 0, j;
  int source, dest;

  int *tempPtr = NULL, *tempArr = (int*)malloc(sizeof(int) * DEFAULT_ARR_SIZE);

  MinHeap* minHeap = createMinHeap(graph->V);

  for (i = 0; i < cases; i++) {
    fscanf(fptr, "%u %u", &source, &dest);

    printf("%d\n", dijkstra(graph, source, dest, distArr, minHeap));
    assert(distArr[dest].steps != INT_MAX);
    assert(distArr[dest].steps > 0);
    if (distArr[dest].steps > (sizeof(tempArr) / sizeof(int))) {
      // resize
      tempPtr = realloc(tempArr, sizeof(int) * distArr[dest].steps);
      assert(tempPtr != NULL);
      tempArr = tempPtr;
      tempPtr = NULL;
    }

    int from;
    for (j = distArr[dest].steps - 1, from = dest; j >= 0; j--) {
      tempArr[j] = from;
      from = distArr[from].from;
    }
    for (j = 0; j < distArr[dest].steps; j++) {
      printf("%d ", tempArr[j]);
    }
    printf("\n");
  }

  destroyMinHeap(minHeap);

  free(tempArr);
  free(tempPtr);

  fclose(fptr);
}

int main(int Argc, char** Argv) {
  // 1. Parse given parameters
  if (Argc < 3) {
    fprintf(stderr, "Must pass two arguments!\n");
    exit(1);
  }
  int numVertex, numEdge;
  Graph* g = read_graph(Argv[1], &numVertex, &numEdge);

  int i = 0;
  DistNode* distArr = (DistNode*)malloc(sizeof(DistNode) * numVertex);
  for (i = 0; i < numVertex; i++) {
    distArr[i].dist = INT_MAX;
    distArr[i].from = 0;
    distArr[i].steps = 0;
  }

  read_query(Argv[2], g, distArr);
  free(distArr);
  destroyGraph(g);

  return 0;
}