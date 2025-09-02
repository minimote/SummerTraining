#include <vector>
#include <queue>
#include <stack>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <memory>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>

namespace AlgorithmLibrary {

// 异常类
class AlgorithmException : public std::runtime_error {
public:
    explicit AlgorithmException(const std::string& message) 
        : std::runtime_error(message) {}
};

// 图基类
class Graph {
protected:
    int vertexCount;
    bool directed;
    
public:
    Graph(int vertices, bool isDirected = false) 
        : vertexCount(vertices), directed(isDirected) {}
    
    virtual ~Graph() = default;
    
    virtual void addEdge(int src, int dest, double weight = 1.0) = 0;
    virtual bool hasEdge(int src, int dest) const = 0;
    virtual double getWeight(int src, int dest) const = 0;
    virtual std::vector<int> getNeighbors(int vertex) const = 0;
    
    int getVertexCount() const { return vertexCount; }
    bool isDirected() const { return directed; }
    
    virtual void print() const = 0;
};

// 邻接矩阵实现的图
class AdjacencyMatrixGraph : public Graph {
private:
    std::vector<std::vector<double>> matrix;
    
public:
    AdjacencyMatrixGraph(int vertices, bool isDirected = false);
    
    void addEdge(int src, int dest, double weight = 1.0) override;
    bool hasEdge(int src, int dest) const override;
    double getWeight(int src, int dest) const override;
    std::vector<int> getNeighbors(int vertex) const override;
    void print() const override;
};

// 邻接表实现的图
class AdjacencyListGraph : public Graph {
private:
    struct Edge {
        int dest;
        double weight;
        Edge(int d, double w) : dest(d), weight(w) {}
    };
    
    std::vector<std::vector<Edge>> adjList;
    
public:
    AdjacencyListGraph(int vertices, bool isDirected = false);
    
    void addEdge(int src, int dest, double weight = 1.0) override;
    bool hasEdge(int src, int dest) const override;
    double getWeight(int src, int dest) const override;
    std::vector<int> getNeighbors(int vertex) const override;
    void print() const override;
};
