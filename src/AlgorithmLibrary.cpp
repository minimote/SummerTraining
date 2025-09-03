#include "../include/AlgorithmLibrary.h"

namespace AlgorithmLibrary
{

    // AdjacencyMatrixGraph 实现
    AdjacencyMatrixGraph::AdjacencyMatrixGraph(int vertices, bool isDirected)
        : Graph(vertices, isDirected), matrix(vertices, std::vector<double>(vertices, AlgorithmUtils::INF))
    {
        for (int i = 0; i < vertices; ++i)
        {
            matrix[i][i] = 0; // 对角线设为0
        }
    }

    void AdjacencyMatrixGraph::addEdge(int src, int dest, double weight)
    {
        if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount)
        {
            throw AlgorithmException("Invalid vertex index");
        }

        matrix[src][dest] = weight;
        if (!directed)
        {
            matrix[dest][src] = weight;
        }
    }

    bool AdjacencyMatrixGraph::hasEdge(int src, int dest) const
    {
        if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount)
        {
            return false;
        }
        return !AlgorithmUtils::approximatelyEqual(matrix[src][dest], AlgorithmUtils::INF);
    }

    double AdjacencyMatrixGraph::getWeight(int src, int dest) const
    {
        if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount)
        {
            throw AlgorithmException("Invalid vertex index");
        }
        return matrix[src][dest];
    }

    std::vector<int> AdjacencyMatrixGraph::getNeighbors(int vertex) const
    {
        if (vertex < 0 || vertex >= vertexCount)
        {
            throw AlgorithmException("Invalid vertex index");
        }

        std::vector<int> neighbors;
        for (int i = 0; i < vertexCount; ++i)
        {
            if (!AlgorithmUtils::approximatelyEqual(matrix[vertex][i], AlgorithmUtils::INF) &&
                !AlgorithmUtils::approximatelyEqual(matrix[vertex][i], 0))
            {
                neighbors.push_back(i);
            }
        }
        return neighbors;
    }

    void AdjacencyMatrixGraph::print() const
    {
        std::cout << "Graph (Adjacency Matrix):" << std::endl;
        for (int i = 0; i < vertexCount; ++i)
        {
            for (int j = 0; j < vertexCount; ++j)
            {
                if (AlgorithmUtils::approximatelyEqual(matrix[i][j], AlgorithmUtils::INF))
                {
                    std::cout << "INF ";
                }
                else
                {
                    std::cout << matrix[i][j] << " ";
                }
            }
            std::cout << std::endl;
        }
    }

    // AdjacencyListGraph 实现
    AdjacencyListGraph::AdjacencyListGraph(int vertices, bool isDirected)
        : Graph(vertices, isDirected), adjList(vertices) {}

    void AdjacencyListGraph::addEdge(int src, int dest, double weight)
    {
        if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount)
        {
            throw AlgorithmException("Invalid vertex index");
        }

        adjList[src].emplace_back(dest, weight);
        if (!directed && src != dest)
        {
            adjList[dest].emplace_back(src, weight);
        }
    }

    bool AdjacencyListGraph::hasEdge(int src, int dest) const
    {
        if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount)
        {
            return false;
        }

        for (const auto &edge : adjList[src])
        {
            if (edge.dest == dest)
            {
                return true;
            }
        }
        return false;
    }

    double AdjacencyListGraph::getWeight(int src, int dest) const
    {
        if (src < 0 || src >= vertexCount || dest < 0 || dest >= vertexCount)
        {
            throw AlgorithmException("Invalid vertex index");
        }

        for (const auto &edge : adjList[src])
        {
            if (edge.dest == dest)
            {
                return edge.weight;
            }
        }
        return AlgorithmUtils::INF;
    }

    std::vector<int> AdjacencyListGraph::getNeighbors(int vertex) const
    {
        if (vertex < 0 || vertex >= vertexCount)
        {
            throw AlgorithmException("Invalid vertex index");
        }

        std::vector<int> neighbors;
        for (const auto &edge : adjList[vertex])
        {
            neighbors.push_back(edge.dest);
        }
        return neighbors;
    }

    void AdjacencyListGraph::print() const
    {
        std::cout << "Graph (Adjacency List):" << std::endl;
        for (int i = 0; i < vertexCount; ++i)
        {
            std::cout << i << ": ";
            for (const auto &edge : adjList[i])
            {
                std::cout << "-> " << edge.dest << "(" << edge.weight << ") ";
            }
            std::cout << std::endl;
        }
    }

    // Graph类中的连通性检测方法实现
    bool Graph::isConnected() const
    {
        // 对于空图或只有一个顶点的图，认为是连通的
        if (getVertexCount() <= 1)
        {
            return true;
        }

        if (isDirected())
        {
            // 对于有向图，检查强连通性
            // 从每个顶点开始都应该能到达所有其他顶点
            for (int start = 0; start < getVertexCount(); ++start)
            {
                std::vector<bool> visited(getVertexCount(), false);
                std::queue<int> queue;

                queue.push(start);
                visited[start] = true;
                int visitedCount = 1;

                while (!queue.empty())
                {
                    int current = queue.front();
                    queue.pop();

                    // 获取当前顶点的所有邻居（出边指向的顶点）
                    std::vector<int> neighbors = getNeighbors(current);
                    for (int neighbor : neighbors)
                    {
                        if (!visited[neighbor])
                        {
                            visited[neighbor] = true;
                            visitedCount++;
                            queue.push(neighbor);
                        }
                    }
                }

                // 如果从这个顶点无法到达所有其他顶点，则不是强连通
                if (visitedCount != getVertexCount())
                {
                    return false;
                }
            }
            return true;
        }
        else
        {
            // 对于无向图，只需从一个顶点开始检查
            std::vector<bool> visited(getVertexCount(), false);
            std::queue<int> queue;

            // 从顶点0开始
            queue.push(0);
            visited[0] = true;
            int visitedCount = 1;

            while (!queue.empty())
            {
                int current = queue.front();
                queue.pop();

                // 获取当前顶点的所有邻居
                std::vector<int> neighbors = getNeighbors(current);
                for (int neighbor : neighbors)
                {
                    if (!visited[neighbor])
                    {
                        visited[neighbor] = true;
                        visitedCount++;
                        queue.push(neighbor);

                        // 如果已经访问了所有顶点，提前结束
                        if (visitedCount == getVertexCount())
                        {
                            return true;
                        }
                    }
                }
            }

            // 如果访问的顶点数等于总顶点数，则图连通
            return visitedCount == getVertexCount();
        }
    }
}
constexpr double AlgorithmLibrary::AlgorithmUtils::INF;
}
