#include "../include/AlgorithmLibrary.h"

namespace AlgorithmLibrary
{

    // AdjacencyMatrixGraph ʵ��
    AdjacencyMatrixGraph::AdjacencyMatrixGraph(int vertices, bool isDirected)
        : Graph(vertices, isDirected), matrix(vertices, std::vector<double>(vertices, AlgorithmUtils::INF))
    {
        for (int i = 0; i < vertices; ++i)
        {
            matrix[i][i] = 0; // �Խ�����Ϊ0
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

    // AdjacencyListGraph ʵ��
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

    // ConnectivityAlgorithms ʵ��
    // ��ͨ�Ե�ʵ��
    bool ConnectivityAlgorithms::isConnected(const Graph &graph)
    {
        if (graph.getVertexCount() == 0)
            return true;

        std::vector<bool> visited(graph.getVertexCount(), false);
        dfsForConnectivity(graph, 0, visited);

        for (bool v : visited)
        {
            if (!v)
                return false;
        }
        return true;
    }

    void ConnectivityAlgorithms::dfsForConnectivity(const Graph &graph, int vertex,
                                                    std::vector<bool> &visited)
    {
        visited[vertex] = true;
        for (int neighbor : graph.getNeighbors(vertex))
        {
            if (!visited[neighbor])
            {
                dfsForConnectivity(graph, neighbor, visited);
            }
        }
    }

    // PathAlgorithms ʵ��

    // ��Դ���·��ʵ��(anjunjie)
    std::vector<std::vector<double>> PathAlgorithms::floydWarshall(const Graph &graph)
    {
        int n = graph.getVertexCount();
        std::vector<std::vector<double>> dist(n, std::vector<double>(n, AlgorithmUtils::INF));

        // ��ʼ���������
        for (int i = 0; i < n; ++i)
        {
            dist[i][i] = 0;
            for (int j : graph.getNeighbors(i))
            {
                dist[i][j] = graph.getWeight(i, j);
            }
        }

        // Floyd-Warshall�㷨����
        for (int k = 0; k < n; ++k)
        {
            for (int i = 0; i < n; ++i)
            {
                if (AlgorithmUtils::approximatelyEqual(dist[i][k], AlgorithmUtils::INF))
                {
                    continue;
                }

                for (int j = 0; j < n; ++j)
                {
                    if (AlgorithmUtils::definitelyLessThan(dist[i][k] + dist[k][j], dist[i][j]))
                    {
                        dist[i][j] = dist[i][k] + dist[k][j];
                    }
                }
            }
        }

        return dist;
    }

    // ��������㷨�͹�������㷨��ʵ��

    std::vector<int> PathAlgorithms::dfsPath(const Graph &graph, int start, int goal)
    {
        int n = graph.getVertexCount();
        if (start < 0 || start >= n || goal < 0 || goal >= n)
        {
            throw AlgorithmException("Invalid start or goal vertex");
        }

        std::vector<bool> visited(n, false);
        std::vector<int> parent(n, -1);
        std::stack<int> stack;

        stack.push(start);
        visited[start] = true;

        while (!stack.empty())
        {
            int current = stack.top();
            stack.pop();

            if (current == goal)
            {
                return reconstructPath(parent, goal);
            }

            for (int neighbor : graph.getNeighbors(current))
            {
                if (!visited[neighbor])
                {
                    visited[neighbor] = true;
                    parent[neighbor] = current;
                    stack.push(neighbor);
                }
            }
        }

        return {}; // û���ҵ�·��
    }

    std::vector<int> PathAlgorithms::bfsPath(const Graph &graph, int start, int goal)
    {
        int n = graph.getVertexCount();
        if (start < 0 || start >= n || goal < 0 || goal >= n)
        {
            throw AlgorithmException("Invalid start or goal vertex");
        }

        std::vector<bool> visited(n, false);
        std::vector<int> parent(n, -1);
        std::queue<int> queue;

        queue.push(start);
        visited[start] = true;

        while (!queue.empty())
        {
            int current = queue.front();
            queue.pop();

            if (current == goal)
            {
                return reconstructPath(parent, goal);
            }

            for (int neighbor : graph.getNeighbors(current))
            {
                if (!visited[neighbor])
                {
                    visited[neighbor] = true;
                    parent[neighbor] = current;
                    queue.push(neighbor);
                }
            }
        }

        return {}; // û���ҵ�·��
    }
    std::vector<int> PathAlgorithms::reconstructPath(const std::vector<int> &prev, int goal)
    {
        std::vector<int> path;
        int current = goal;

        while (current != -1)
        {
            path.push_back(current);
            current = prev[current];
        }

        std::reverse(path.begin(), path.end());
        return path;
    }

    // AlgorithmUtils ʵ��
    bool AlgorithmUtils::approximatelyEqual(double a, double b, double epsilon)
    {
        return std::abs(a - b) <= epsilon;
    }

    bool AlgorithmUtils::definitelyLessThan(double a, double b, double epsilon)
    {
        return (b - a) > epsilon;
    }

    bool AlgorithmUtils::definitelyGreaterThan(double a, double b, double epsilon)
    {
        return (a - b) > epsilon;
    }
    constexpr double AlgorithmLibrary::AlgorithmUtils::INF;

    // EulerianAlgorithms ʵ��
    bool EulerianAlgorithms::hasEulerianCircuit(const Graph &graph)
    {
        int n = graph.getVertexCount();
        if (n == 0)
            return true;

        // ���ͼ�Ƿ���ͨ
        if (!ConnectivityAlgorithms::isConnected(graph))
        {
            return false;
        }

        if (graph.isDirected())
        {
            // ����ͼ���
            for (int i = 0; i < n; ++i)
            {
                auto degrees = getInDegreeAndOutDegree(graph, i);
                if (degrees.first != degrees.second)
                {
                    return false;
                }
            }
        }
        else
        {
            // ����ͼ���
            for (int i = 0; i < n; ++i)
            {
                if (getDegree(graph, i) % 2 != 0)
                {
                    return false;
                }
            }
        }

        return true;
    }

    bool EulerianAlgorithms::hasEulerianPath(const Graph &graph)
    {
        int n = graph.getVertexCount();
        if (n == 0)
            return true;

        // ���ͼ�Ƿ���ͨ
        if (!ConnectivityAlgorithms::isConnected(graph))
        {
            return false;
        }

        if (graph.isDirected())
        {
            // ����ͼ���
            int startVertices = 0; // ���ȱ���ȴ�1�Ķ�����
            int endVertices = 0;   // ��ȱȳ��ȴ�1�Ķ�����

            for (int i = 0; i < n; ++i)
            {
                auto degrees = getInDegreeAndOutDegree(graph, i);
                int inDegree = degrees.first;
                int outDegree = degrees.second;

                if (outDegree - inDegree == 1)
                {
                    startVertices++;
                }
                else if (inDegree - outDegree == 1)
                {
                    endVertices++;
                }
                else if (inDegree != outDegree)
                {
                    return false;
                }
            }

            return (startVertices == 0 && endVertices == 0) ||
                   (startVertices == 1 && endVertices == 1);
        }
        else
        {
            // ����ͼ���
            int oddDegreeVertices = 0;
            for (int i = 0; i < n; ++i)
            {
                if (getDegree(graph, i) % 2 != 0)
                {
                    oddDegreeVertices++;
                }
            }

            return oddDegreeVertices == 0 || oddDegreeVertices == 2;
        }
    }

    std::vector<int> EulerianAlgorithms::findEulerianCircuit(const Graph &graph)
    {
        if (!hasEulerianCircuit(graph))
        {
            throw AlgorithmException("Graph does not have an Eulerian circuit");
        }

        int n = graph.getVertexCount();
        if (n == 0)
            return {};

        std::vector<int> circuit;
        std::vector<int> currentPath;
        std::vector<std::vector<bool>> visitedEdges(n, std::vector<bool>(n, false));

        // �Ӷ���0��ʼ����
        dfsEulerian(graph, 0, visitedEdges, circuit, currentPath);

        return circuit;
    }

    std::vector<int> EulerianAlgorithms::findEulerianPath(const Graph &graph)
    {
        if (!hasEulerianPath(graph))
        {
            throw AlgorithmException("Graph does not have an Eulerian path");
        }

        int n = graph.getVertexCount();
        if (n == 0)
            return {};

        std::vector<int> path;
        std::vector<int> currentPath;
        std::vector<std::vector<bool>> visitedEdges(n, std::vector<bool>(n, false));

        int startVertex = 0;
        if (graph.isDirected())
        {
            // Ѱ����ʼ���㣨���ȴ�����ȵĶ��㣩
            for (int i = 0; i < n; ++i)
            {
                auto degrees = getInDegreeAndOutDegree(graph, i);
                if (degrees.second > degrees.first)
                {
                    startVertex = i;
                    break;
                }
            }
        }
        else
        {
            // Ѱ����ȶ�����Ϊ���
            for (int i = 0; i < n; ++i)
            {
                if (getDegree(graph, i) % 2 != 0)
                {
                    startVertex = i;
                    break;
                }
            }
        }

        dfsEulerian(graph, startVertex, visitedEdges, path, currentPath);

        return path;
    }

    void EulerianAlgorithms::dfsEulerian(const Graph &graph, int vertex, std::vector<std::vector<bool>> &visitedEdges,
                                         std::vector<int> &circuit, std::vector<int> &currentPath)
    {
        // ������ǰ����������ھ�
        std::vector<int> neighbors = graph.getNeighbors(vertex);

        for (int neighbor : neighbors)
        {
            // �����Ƿ��Ѿ�������
            if (!visitedEdges[vertex][neighbor])
            {
                // ��Ǳ�Ϊ�ѷ���
                visitedEdges[vertex][neighbor] = true;
                // ��������ͼ��Ҫͬʱ��Ƿ����
                if (!graph.isDirected())
                {
                    visitedEdges[neighbor][vertex] = true;
                }

                // �ݹ����
                dfsEulerian(graph, neighbor, visitedEdges, circuit, currentPath);
            }
        }

        // û�и���δ���ʵıߣ���������뵽·����
        circuit.push_back(vertex);
    }

    int EulerianAlgorithms::getDegree(const Graph &graph, int vertex)
    {
        return static_cast<int>(graph.getNeighbors(vertex).size());
    }

    std::pair<int, int> EulerianAlgorithms::getInDegreeAndOutDegree(const Graph &graph, int vertex)
    {
        // ����ֵ: (���, ����)
        int inDegree = 0;
        int outDegree = 0;

        outDegree = getDegree(graph, vertex);

        // �������
        for (int i = 0; i < graph.getVertexCount(); ++i)
        {
            if (i != vertex)
            {
                std::vector<int> neighbors = graph.getNeighbors(i);
                if (std::find(neighbors.begin(), neighbors.end(), vertex) != neighbors.end())
                {
                    inDegree++;
                }
            }
        }

        return std::make_pair(inDegree, outDegree);
    }

} // namespace AlgorithmLibrary
