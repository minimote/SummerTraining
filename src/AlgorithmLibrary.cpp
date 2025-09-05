#include "../include/AlgorithmLibrary.h"

namespace AlgorithmLibrary {

// AdjacencyMatrixGraph 实现-安钧杰 
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

// AdjacencyListGraph 实现-安钧杰 
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

// ConnectivityAlgorithms 实现
 
bool ConnectivityAlgorithms::isConnected(const Graph& graph) {
    if (graph.getVertexCount() == 0) return true;
    
    std::vector<bool> visited(graph.getVertexCount(), false);
    dfsForConnectivity(graph, 0, visited);
    
    for (bool v : visited) {
        if (!v) return false;
    }
    return true;
}

std::vector<std::vector<int>> ConnectivityAlgorithms::findConnectedComponents(const Graph& graph) {
    std::vector<std::vector<int>> components;
    std::vector<bool> visited(graph.getVertexCount(), false);
    
    for (int i = 0; i < graph.getVertexCount(); ++i) {
        if (!visited[i]) {
            std::vector<int> component;
            std::stack<int> stack;
            stack.push(i);
            visited[i] = true;
            
            while (!stack.empty()) {
                int vertex = stack.top();
                stack.pop();
                component.push_back(vertex);
                
                for (int neighbor : graph.getNeighbors(vertex)) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        stack.push(neighbor);
                    }
                }
            }
            
            components.push_back(component);
        }
    }
    
    return components;
}

std::vector<int> ConnectivityAlgorithms::findArticulationPoints(const Graph& graph) {
    int n = graph.getVertexCount();
    std::vector<int> disc(n, -1);
    std::vector<int> low(n, -1);
    std::vector<int> parent(n, -1);
    std::vector<bool> ap(n, false);
    int time = 0;
    
    for (int i = 0; i < n; ++i) {
        if (disc[i] == -1) {
            dfsForArticulationPoints(graph, i, disc, low, parent, ap, time);
        }
    }
    
    std::vector<int> result;
    for (int i = 0; i < n; ++i) {
        if (ap[i]) {
            result.push_back(i);
        }
    }
    
    return result;
}

std::vector<std::pair<int, int>> ConnectivityAlgorithms::findBridges(const Graph& graph) {
    int n = graph.getVertexCount();
    std::vector<int> disc(n, -1);
    std::vector<int> low(n, -1);
    std::vector<int> parent(n, -1);
    std::vector<std::pair<int, int>> bridges;
    int time = 0;
    
    for (int i = 0; i < n; ++i) {
        if (disc[i] == -1) {
            dfsForBridges(graph, i, disc, low, parent, bridges, time);
        }
    }
    
    return bridges;
}

void ConnectivityAlgorithms::dfsForConnectivity(const Graph& graph, int vertex, 
                                               std::vector<bool>& visited) {
    visited[vertex] = true;
    for (int neighbor : graph.getNeighbors(vertex)) {
        if (!visited[neighbor]) {
            dfsForConnectivity(graph, neighbor, visited);
        }
    }
}

void ConnectivityAlgorithms::dfsForArticulationPoints(const Graph& graph, int vertex, 
                                                    std::vector<int>& disc, std::vector<int>& low, 
                                                    std::vector<int>& parent, std::vector<bool>& ap, 
                                                    int& time) {
    disc[vertex] = low[vertex] = ++time;
    int children = 0;
    
    for (int neighbor : graph.getNeighbors(vertex)) {
        if (disc[neighbor] == -1) {
            children++;
            parent[neighbor] = vertex;
            dfsForArticulationPoints(graph, neighbor, disc, low, parent, ap, time);
            
            low[vertex] = std::min(low[vertex], low[neighbor]);
            
            // 根节点且有多个子节点
            if (parent[vertex] == -1 && children > 1) {
                ap[vertex] = true;
            }
            
            // 非根节点且low[neighbor] >= disc[vertex]
            if (parent[vertex] != -1 && low[neighbor] >= disc[vertex]) {
                ap[vertex] = true;
            }
        } else if (neighbor != parent[vertex]) {
            low[vertex] = std::min(low[vertex], disc[neighbor]);
        }
    }
}

void ConnectivityAlgorithms::dfsForBridges(const Graph& graph, int vertex, 
                                         std::vector<int>& disc, std::vector<int>& low, 
                                         std::vector<int>& parent, 
                                         std::vector<std::pair<int, int>>& bridges, 
                                         int& time) {
    disc[vertex] = low[vertex] = ++time;
    
    for (int neighbor : graph.getNeighbors(vertex)) {
        if (disc[neighbor] == -1) {
            parent[neighbor] = vertex;
            dfsForBridges(graph, neighbor, disc, low, parent, bridges, time);
            
            low[vertex] = std::min(low[vertex], low[neighbor]);
            
            if (low[neighbor] > disc[vertex]) {
                bridges.emplace_back(vertex, neighbor);
            }
        } else if (neighbor != parent[vertex]) {
            low[vertex] = std::min(low[vertex], disc[neighbor]);
        }
    }
}

// PathAlgorithms 实现
std::vector<double> PathAlgorithms::dijkstra(const Graph& graph, int source) {
    int n = graph.getVertexCount();
    if (source < 0 || source >= n) {
        throw AlgorithmException("Invalid source vertex");
    }
    
    std::vector<double> dist(n, AlgorithmUtils::INF);
    dist[source] = 0;
    
    // 使用优先队列（最小堆）
    using Pair = std::pair<double, int>;
    std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair>> pq;
    pq.emplace(0, source);
    
    while (!pq.empty()) {
        double d = pq.top().first;
        int u = pq.top().second;
        pq.pop();
        
        // 如果找到更短的路径已经处理过，跳过
        if (d > dist[u]) {
            continue;
        }
        
        for (int v : graph.getNeighbors(u)) {
            double weight = graph.getWeight(u, v);
            if (AlgorithmUtils::definitelyGreaterThan(dist[u] + weight, dist[v])) {
                continue;
            }
            
            if (AlgorithmUtils::definitelyLessThan(dist[u] + weight, dist[v])) {
                dist[v] = dist[u] + weight;
                pq.emplace(dist[v], v);
            }
        }
    }
    
    return dist;
}

std::vector<double> PathAlgorithms::bellmanFord(const Graph& graph, int source) {
    int n = graph.getVertexCount();
    if (source < 0 || source >= n) {
        throw AlgorithmException("Invalid source vertex");
    }
    
    std::vector<double> dist(n, AlgorithmUtils::INF);
    dist[source] = 0;
    
    // 松弛操作执行V-1次
    for (int i = 1; i < n; ++i) {
        for (int u = 0; u < n; ++u) {
            if (AlgorithmUtils::approximatelyEqual(dist[u], AlgorithmUtils::INF)) {
                continue;
            }
            
            for (int v : graph.getNeighbors(u)) {
                double weight = graph.getWeight(u, v);
                if (AlgorithmUtils::definitelyLessThan(dist[u] + weight, dist[v])) {
                    dist[v] = dist[u] + weight;
                }
            }
        }
    }
    
    // 检查负权环
    for (int u = 0; u < n; ++u) {
        if (AlgorithmUtils::approximatelyEqual(dist[u], AlgorithmUtils::INF)) {
            continue;
        }
        
        for (int v : graph.getNeighbors(u)) {
            double weight = graph.getWeight(u, v);
            if (AlgorithmUtils::definitelyLessThan(dist[u] + weight, dist[v])) {
                throw AlgorithmException("Graph contains a negative-weight cycle");
            }
        }
    }
    
    return dist;
}

std::vector<std::vector<double>> PathAlgorithms::floydWarshall(const Graph& graph) {
    int n = graph.getVertexCount();
    std::vector<std::vector<double>> dist(n, std::vector<double>(n, AlgorithmUtils::INF));
    
    // 初始化距离矩阵
    for (int i = 0; i < n; ++i) {
        dist[i][i] = 0;
        for (int j : graph.getNeighbors(i)) {
            dist[i][j] = graph.getWeight(i, j);
        }
    }
    
    // Floyd-Warshall算法核心
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            if (AlgorithmUtils::approximatelyEqual(dist[i][k], AlgorithmUtils::INF)) {
                continue;
            }
            
            for (int j = 0; j < n; ++j) {
                if (AlgorithmUtils::definitelyLessThan(dist[i][k] + dist[k][j], dist[i][j])) {
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }
    
    return dist;
}

std::vector<int> PathAlgorithms::aStar(const Graph& graph, int start, int goal,
                                      double (*heuristic)(int, int)) {
    int n = graph.getVertexCount();
    if (start < 0 || start >= n || goal < 0 || goal >= n) {
        throw AlgorithmException("Invalid start or goal vertex");
    }
    
    std::vector<double> gScore(n, AlgorithmUtils::INF);
    std::vector<double> fScore(n, AlgorithmUtils::INF);
    std::vector<int> cameFrom(n, -1);
    
    gScore[start] = 0;
    fScore[start] = heuristic(start, goal);
    
    using Pair = std::pair<double, int>;
    std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair>> openSet;
    openSet.emplace(fScore[start], start);
    
    while (!openSet.empty()) {
        int current = openSet.top().second;
        openSet.pop();
        
        if (current == goal) {
            return reconstructPath(cameFrom, goal);
        }
        
        for (int neighbor : graph.getNeighbors(current)) {
            double tentativeGScore = gScore[current] + graph.getWeight(current, neighbor);
            
            if (AlgorithmUtils::definitelyLessThan(tentativeGScore, gScore[neighbor])) {
                cameFrom[neighbor] = current;
                gScore[neighbor] = tentativeGScore;
                fScore[neighbor] = gScore[neighbor] + heuristic(neighbor, goal);
                openSet.emplace(fScore[neighbor], neighbor);
            }
        }
    }
    
    return {}; // 没有找到路径
}

std::vector<int> PathAlgorithms::dfsPath(const Graph& graph, int start, int goal) {
    int n = graph.getVertexCount();
    if (start < 0 || start >= n || goal < 0 || goal >= n) {
        throw AlgorithmException("Invalid start or goal vertex");
    }
    
    std::vector<bool> visited(n, false);
    std::vector<int> parent(n, -1);
    std::stack<int> stack;
    
    stack.push(start);
    visited[start] = true;
    
    while (!stack.empty()) {
        int current = stack.top();
        stack.pop();
        
        if (current == goal) {
            return reconstructPath(parent, goal);
        }
        
        for (int neighbor : graph.getNeighbors(current)) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                parent[neighbor] = current;
                stack.push(neighbor);
            }
        }
    }
    
    return {}; // 没有找到路径
}

std::vector<int> PathAlgorithms::bfsPath(const Graph& graph, int start, int goal) {
    int n = graph.getVertexCount();
    if (start < 0 || start >= n || goal < 0 || goal >= n) {
        throw AlgorithmException("Invalid start or goal vertex");
    }
    
    std::vector<bool> visited(n, false);
    std::vector<int> parent(n, -1);
    std::queue<int> queue;
    
    queue.push(start);
    visited[start] = true;
    
    while (!queue.empty()) {
        int current = queue.front();
        queue.pop();
        
        if (current == goal) {
            return reconstructPath(parent, goal);
        }
        
        for (int neighbor : graph.getNeighbors(current)) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                parent[neighbor] = current;
                queue.push(neighbor);
            }
        }
    }
    
    return {}; // 没有找到路径
}

std::vector<int> PathAlgorithms::reconstructPath(const std::vector<int>& prev, int goal) {
    std::vector<int> path;
    int current = goal;
    
    while (current != -1) {
        path.push_back(current);
        current = prev[current];
    }
    
    std::reverse(path.begin(), path.end());
    return path;
}

// AlgorithmUtils 实现
bool AlgorithmUtils::approximatelyEqual(double a, double b, double epsilon) {
    return std::abs(a - b) <= epsilon;
}

bool AlgorithmUtils::definitelyLessThan(double a, double b, double epsilon) {
    return (b - a) > epsilon;
}

bool AlgorithmUtils::definitelyGreaterThan(double a, double b, double epsilon) {
    return (a - b) > epsilon;
}
constexpr double AlgorithmLibrary::AlgorithmUtils::INF;
} // namespace AlgorithmLibrary
