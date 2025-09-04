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

// �쳣��
class AlgorithmException : public std::runtime_error {
public:
    explicit AlgorithmException(const std::string& message) 
        : std::runtime_error(message) {}
};

// ͼ����
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

// �ڽӾ���ʵ�ֵ�ͼ
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

// �ڽӱ�ʵ�ֵ�ͼ
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

// ��ͨ���㷨��
class ConnectivityAlgorithms {
public:
    // ���ͼ�Ƿ���ͨ
    static bool isConnected(const Graph& graph);
    
    // ������ͨ����
    static std::vector<std::vector<int>> findConnectedComponents(const Graph& graph);
    
    // ���ҹؽڵ㣨��㣩
    static std::vector<int> findArticulationPoints(const Graph& graph);
    
    // ������
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


// ŷ���㷨
class EulerianAlgorithms {
public:
    // �ж�ͼ�Ƿ���ŷ����·
    static bool hasEulerianCircuit(const Graph& graph);
    
    // �ж�ͼ�Ƿ���ŷ��·��
    static bool hasEulerianPath(const Graph& graph);
    
    // ����ŷ����·
    static std::vector<int> findEulerianCircuit(const Graph& graph);
    
    // ����ŷ��·��
    static std::vector<int> findEulerianPath(const Graph& graph);

private:
    // ������ȱ�����������
    static void dfsEulerian(const Graph& graph, int vertex, std::vector<std::vector<bool>>& visitedEdges, 
                           std::vector<int>& circuit, std::vector<int>& currentPath);
    
    // ��ȡ�������
    static int getDegree(const Graph& graph, int vertex);
    
    // ��ȡ������Ⱥͳ���(����ͼ)
    static std::pair<int, int> getInDegreeAndOutDegree(const Graph& graph, int vertex);
};


// ·���㷨��
class PathAlgorithms {
public:
    // Dijkstra�㷨 - ��Դ���·��
    static std::vector<double> dijkstra(const Graph& graph, int source);
    
    // Bellman-Ford�㷨 - ����Ȩ��
    static std::vector<double> bellmanFord(const Graph& graph, int source);
    
    // Floyd-Warshall�㷨 - ���нڵ�����·��
    static std::vector<std::vector<double>> floydWarshall(const Graph& graph);
    
    // A*�㷨 - ����ʽ����
    static std::vector<int> aStar(const Graph& graph, int start, int goal,
                                 double (*heuristic)(int, int));
    
    // �����������·��
    static std::vector<int> dfsPath(const Graph& graph, int start, int goal);
    
    // �����������·��
    static std::vector<int> bfsPath(const Graph& graph, int start, int goal);
    
    // �ؽ�·���ĸ�������
    static std::vector<int> reconstructPath(const std::vector<int>& prev, int goal);
};

// ���ߺ���
class AlgorithmUtils {
public:
    static constexpr double INF = std::numeric_limits<double>::infinity();
    
    static bool approximatelyEqual(double a, double b, double epsilon = 1e-9);
    static bool definitelyLessThan(double a, double b, double epsilon = 1e-9);
    static bool definitelyGreaterThan(double a, double b, double epsilon = 1e-9);
};

} // namespace AlgorithmLibrary

#endif // ALGORITHM_LIBRARY_H
