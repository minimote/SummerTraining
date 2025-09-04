#ifndef ALGORITHM_LIBRARY_H
#define ALGORITHM_LIBRARY_H

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

// 连通性算法类
class ConnectivityAlgorithms {
public:
    // 检查图是否连通
    static bool isConnected(const Graph& graph);
    
    // 查找连通分量
    static std::vector<std::vector<int>> findConnectedComponents(const Graph& graph);
    
    // 查找关节点（割点）
    static std::vector<int> findArticulationPoints(const Graph& graph);
    
    // 查找桥
    static std::vector<std::pair<int, int>> findBridges(const Graph& graph);
    
private:
    static void dfsForConnectivity(const Graph& graph, int vertex, 
                                  std::vector<bool>& visited);
    static void dfsForArticulationPoints(const Graph& graph, int vertex, 
                                       std::vector<int>& disc, std::vector<int>& low, 
                                       std::vector<int>& parent, std::vector<bool>& ap, 
                                       int& time);
    static void dfsForBridges(const Graph& graph, int vertex, 
                            std::vector<int>& disc, std::vector<int>& low, 
                            std::vector<int>& parent, 
                            std::vector<std::pair<int, int>>& bridges, int& time);
};


// 欧拉算法
class EulerianAlgorithms {
public:
    // 判断图是否有欧拉回路
    static bool hasEulerianCircuit(const Graph& graph);
    
    // 判断图是否有欧拉路径
    static bool hasEulerianPath(const Graph& graph);
    
    // 查找欧拉回路
    static std::vector<int> findEulerianCircuit(const Graph& graph);
    
    // 查找欧拉路径
    static std::vector<int> findEulerianPath(const Graph& graph);

private:
    // 深度优先遍历辅助函数
    static void dfsEulerian(const Graph& graph, int vertex, std::vector<std::vector<bool>>& visitedEdges, 
                           std::vector<int>& circuit, std::vector<int>& currentPath);
    
    // 获取顶点度数
    static int getDegree(const Graph& graph, int vertex);
    
    // 获取顶点入度和出度(有向图)
    static std::pair<int, int> getInDegreeAndOutDegree(const Graph& graph, int vertex);
};


// 路径算法类
class PathAlgorithms {
public:
    // Dijkstra算法 - 单源最短路径
    static std::vector<double> dijkstra(const Graph& graph, int source);
    
    // Bellman-Ford算法 - 处理负权边
    static std::vector<double> bellmanFord(const Graph& graph, int source);
    
    // Floyd-Warshall算法 - 所有节点对最短路径
    static std::vector<std::vector<double>> floydWarshall(const Graph& graph);
    
    // A*算法 - 启发式搜索
    static std::vector<int> aStar(const Graph& graph, int start, int goal,
                                 double (*heuristic)(int, int));
    
    // 深度优先搜索路径
    static std::vector<int> dfsPath(const Graph& graph, int start, int goal);
    
    // 广度优先搜索路径
    static std::vector<int> bfsPath(const Graph& graph, int start, int goal);
    
    // 重建路径的辅助函数
    static std::vector<int> reconstructPath(const std::vector<int>& prev, int goal);
};

// 工具函数
class AlgorithmUtils {
public:
    static constexpr double INF = std::numeric_limits<double>::infinity();
    
    static bool approximatelyEqual(double a, double b, double epsilon = 1e-9);
    static bool definitelyLessThan(double a, double b, double epsilon = 1e-9);
    static bool definitelyGreaterThan(double a, double b, double epsilon = 1e-9);
};

} // namespace AlgorithmLibrary

#endif // ALGORITHM_LIBRARY_H
