#include "../include/AlgorithmLibrary.h"
#include <iostream>

using namespace AlgorithmLibrary;

// A*算法的启发式函数示例（曼哈顿距离）
double manhattanHeuristic(int a, int b) {
    // 这里假设顶点编号对应坐标，实际应用中需要根据具体问题设计
    return std::abs(a - b);
}

int main() {
    try {
        // 创建图
        AdjacencyListGraph graph(6, false);
        graph.addEdge(0, 1, 7);
        graph.addEdge(0, 2, 9);
        graph.addEdge(0, 5, 14);
        graph.addEdge(1, 2, 10);
        graph.addEdge(1, 3, 15);
        graph.addEdge(2, 3, 11);
        graph.addEdge(2, 5, 2);
        graph.addEdge(3, 4, 6);
        graph.addEdge(4, 5, 9);
        
        graph.print();
        
        // 连通性检查
        std::cout << "Graph is connected: " << ConnectivityAlgorithms::isConnected(graph) << std::endl;
        
        // 查找连通分量
        auto components = ConnectivityAlgorithms::findConnectedComponents(graph);
        std::cout << "Connected components: " << components.size() << std::endl;
        for (size_t i = 0; i < components.size(); ++i) {
            std::cout << "Component " << i << ": ";
            for (int vertex : components[i]) {
                std::cout << vertex << " ";
            }
            std::cout << std::endl;
        }
        
        // 最短路径算法
        std::cout << "\nDijkstra from vertex 0:" << std::endl;
        auto distances = PathAlgorithms::dijkstra(graph, 0);
        for (size_t i = 0; i < distances.size(); ++i) {
            std::cout << "Distance to " << i << ": " << distances[i] << std::endl;
        }
        
        // BFS路径查找
        std::cout << "\nBFS path from 0 to 4:" << std::endl;
        auto path = PathAlgorithms::bfsPath(graph, 0, 4);
        for (int vertex : path) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
        
        // A*算法
        std::cout << "\nA* path from 0 to 4:" << std::endl;
        auto aStarPath = PathAlgorithms::aStar(graph, 0, 4, manhattanHeuristic);
        for (int vertex : aStarPath) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
        
    } catch (const AlgorithmException& e) {
        std::cerr << "Algorithm error: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    
    return 0;
}
