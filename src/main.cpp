#include "../include/AlgorithmLibrary.h"
#include <iostream>

using namespace AlgorithmLibrary;
using namespace std;
// A*算法的启发式函数示例（曼哈顿距离）
double manhattanHeuristic(int a, int b) {
    // 这里假设顶点编号对应坐标，实际应用中需要根据具体问题设计
    return std::abs(a - b);
}

int main() {
    
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
        std::cout << "Graph is connected: " << ConnectivityAlgorithms::isConnected(graph) << std::endl;
        
        
    
    return 0;

}
