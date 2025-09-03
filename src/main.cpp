#include "../include/AlgorithmLibrary.h"
#include <iostream>

using namespace AlgorithmLibrary;
using namespace std;
// A*�㷨������ʽ����ʾ���������پ��룩
double manhattanHeuristic(int a, int b) {
    // ������趥���Ŷ�Ӧ���꣬ʵ��Ӧ������Ҫ���ݾ����������
    return std::abs(a - b);
}

int main() {
    
        // ����ͼ
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
