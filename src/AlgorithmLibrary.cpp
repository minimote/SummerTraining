#include "../include/AlgorithmLibrary.h"

namespace AlgorithmLibrary {

// AdjacencyMatrixGraph 实现
AdjacencyMatrixGraph::AdjacencyMatrixGraph(int vertices, bool isDirected)
    : Graph(vertices, isDirected), matrix(vertices, std::vector<double>(vertices, AlgorithmUtils::INF)) {
    for (int i = 0; i < vertices; ++i) {
        matrix[i][i] = 0; // 对角线设为0
    }
}

void AdjacencyMatrixGraph::addEdge(int src, int dest, double weight) {
    if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount) {
        throw AlgorithmException("Invalid vertex index");
    }
    
    matrix[src][dest] = weight;
    if (!directed) {
        matrix[dest][src] = weight;
    }
}

bool AdjacencyMatrixGraph::hasEdge(int src, int dest) const {
    if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount) {
        return false;
    }
    return !AlgorithmUtils::approximatelyEqual(matrix[src][dest], AlgorithmUtils::INF);
}

double AdjacencyMatrixGraph::getWeight(int src, int dest) const {
    if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount) {
        throw AlgorithmException("Invalid vertex index");
    }
    return matrix[src][dest];
}

std::vector<int> AdjacencyMatrixGraph::getNeighbors(int vertex) const {
    if (vertex < 0 || vertex >= vertexCount) {
        throw AlgorithmException("Invalid vertex index");
    }
    
    std::vector<int> neighbors;
    for (int i = 0; i < vertexCount; ++i) {
        if (!AlgorithmUtils::approximatelyEqual(matrix[vertex][i], AlgorithmUtils::INF) && 
            !AlgorithmUtils::approximatelyEqual(matrix[vertex][i], 0)) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

void AdjacencyMatrixGraph::print() const {
    std::cout << "Graph (Adjacency Matrix):" << std::endl;
    for (int i = 0; i < vertexCount; ++i) {
        for (int j = 0; j < vertexCount; ++j) {
            if (AlgorithmUtils::approximatelyEqual(matrix[i][j], AlgorithmUtils::INF)) {
                std::cout << "INF ";
            } else {
                std::cout << matrix[i][j] << " ";
            }
        }
        std::cout << std::endl;
    }
}

// AdjacencyListGraph 实现
AdjacencyListGraph::AdjacencyListGraph(int vertices, bool isDirected)
    : Graph(vertices, isDirected), adjList(vertices) {}

void AdjacencyListGraph::addEdge(int src, int dest, double weight) {
    if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount) {
        throw AlgorithmException("Invalid vertex index");
    }
    
    adjList[src].emplace_back(dest, weight);
    if (!directed && src != dest) {
        adjList[dest].emplace_back(src, weight);
    }
}

bool AdjacencyListGraph::hasEdge(int src, int dest) const {
    if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount) {
        return false;
    }
    
    for (const auto& edge : adjList[src]) {
        if (edge.dest == dest) {
            return true;
        }
    }
    return false;
}

double AdjacencyListGraph::getWeight(int src, int dest) const {
    if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount) {
        throw AlgorithmException("Invalid vertex index");
    }
    
    for (const auto& edge : adjList[src]) {
        if (edge.dest == dest) {
            return edge.weight;
        }
    }
    return AlgorithmUtils::INF;
}

std::vector<int> AdjacencyListGraph::getNeighbors(int vertex) const {
    if (vertex < 0 || vertex >= vertexCount) {
        throw AlgorithmException("Invalid vertex index");
    }
    
    std::vector<int> neighbors;
    for (const auto& edge : adjList[vertex]) {
        neighbors.push_back(edge.dest);
    }
    return neighbors;
}

void AdjacencyListGraph::print() const {
    std::cout << "Graph (Adjacency List):" << std::endl;
    for (int i = 0; i < vertexCount; ++i) {
        std::cout << i << ": ";
        for (const auto& edge : adjList[i]) {
            std::cout << "-> " << edge.dest << "(" << edge.weight << ") ";
        }
        std::cout << std::endl;
    }
}
