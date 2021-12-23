#include <bits/stdc++.h>

using namespace std;

/// For disjoint sets
class Forest
{
private:
    int parent, size;
public:
    static vector<Forest> MakeForest(const int& n)
    {
        vector<Forest> F(n + 1);
        for (int i = 1; i <= n; ++i)
        {
            F[i].parent = i;
            F[i].size = 1;
        }

        return F;
    }

    static int FindParent(int x, vector<Forest>& F)
    {
        while (F[x].parent != x)
        {
            F[x].parent = F[ F[x].parent].parent;
            F[x].size = F[F[x].parent].size;
            x = F[x].parent;
        }
        return x;
    }

    static void UnionForest(int x, int y, vector<Forest>& F)
    {
        x = FindParent(x, F);
        y = FindParent(y, F);

        if (x == y) return;

        if (F[x].size  > F[y].size)
        {
            F[y].parent = x;
            F[x].size += F[y].size;
        }
        else
        {
            F[x].parent = y;
            F[x].size += F[x].size;
        }
    }

    static bool CheckForest(int x, int y, vector<Forest>& F)
    {
        if (FindParent(x, F) == FindParent(y, F))
            return 1;
        return 0;
    }
};

class Edge
{
    //class that helps with the definition of edges, with costs and or capacity
public:
    int from, to;
    int weight;
    int capacity;
    Edge(const int& f=0, const int& w=0, const int& wg = 0, const int& cap = 0)
    {
        from = f;
        to = w;
        weight = wg;
        capacity = cap;
    }
     bool operator()(const Edge &left, const Edge& right) // for priority queue compare
    {
        return left.weight > right.weight;
    }
};

/// Function used to sort a vector of pairs descending by the second component
bool SortDescPair(const pair<int,int>& i, const pair<int,int>& j)
{
    return i.second > j.second;
}

class Graph
{
 private:
    int _nr_vertex; // Number of vertexes
    int _nr_edge; // Number of edges
    bool _oriented; // Mark if the graph is oriented
    bool _weighted; // Mark if the graph is weighted
    bool _flow; // Mark if it is a flow network
    vector<vector<Edge>> _adjacency;

    vector<int> _popEdges(stack<Edge>& , const Edge&);
    void _popVertex(const int&, stack<int>&, int&, vector<vector<int>>&, vector<bool>&);
    int  _bfs(const int&, int&);
    vector<int> _dist(const int&);
    void _dfs(const int&, const int&, vector<int>&, stack<int>&);
    void _leveledDfs(const int&, vector<int>& , vector<int>&, vector<int>& , stack<Edge>&, vector<vector<int>>&, vector<vector<int>>&);
    void _Tarjan(const int&, int&, vector<int>&, vector<int>&, stack<int>&, vector<bool>& , int&, vector<vector<int>>&);
    bool _HavelHakimi(const int&, const int&, vector<pair<int,int>>& );
    void _addEdge(const Edge&);
    void _resize(const int&, const int&, const bool&, const bool&, const bool&);
    int _flowBfs(const vector<vector<int>>&, const int&, const int&, vector<int>&, vector<vector<int>>&, vector<vector<int>>&);
    void _initializeAdMat(vector<vector<int>>&);

 public:
    Graph(const bool& oriented = 0, const bool& weighted = 0, const bool &flow = 0);

    void AddEdge(const int& from, const int& to, const int& cost = 0, const int& capacity = 0);
    void clear();
    void resize(const int& new_nr_nodes);
    int nrVertex();
    int nrEdges();
    void MakeOriented();
    void MakeNonOriented();
    void MakeWeighted();
    void MakeNonWeighted();
    void MakeFlowNetwork();
    void MakeNotFlowNetwork();

    friend istream& operator>>(istream& in, Graph& g);
    friend ostream& operator<<(ostream& out, const Graph& g);

    vector<int> DF(const int&);
    vector<int> BF(const int&);
    int NrConnectedComponents();
    vector<int> DistanceFrom(const int&);
    pair<int,vector<Edge>> Prim(const int&);
    vector<vector<int>> RoyFloyd(vector<vector<int>>& cost_matrix);
    vector<int> TopologicalSort();
    pair<int,vector<vector<int>>> BiconnectedComponents();
    vector<vector<int>> CriticalConnections();
    pair<int, vector<vector<int>>> StronglyConnectedComponents();
    void HavelHakimi(istream&);
    static bool VerifyHavelHakimi(const int&, const int&, vector<pair<int,int>>&);
    tuple<int, int, vector<Edge>> MinimumSpanningTree();
    vector<int> Dijkstra(const int&);
    vector<int> BellmanFord(const int&);
    int  MaxFlow(const int&, const int &);
    int Diameter();
    vector<int> EulerianCycle();
};

istream& operator>>(istream& in, Graph& g)
{
    in >> g._nr_vertex >> g._nr_edge;

    g._adjacency.resize(g._nr_vertex + 1);
    int x, y, c = 0, capacity = 0;
    for(int i = 0; i < g._nr_edge; ++i)
    {
        in >> x >> y;
        if (g._weighted == 1)
        {
            in >> c;
        }
        if (g._flow == 1)
        {
            in >> capacity;
        }
        g._addEdge({x, y, c, capacity});
    }

    return in;
}

ostream& operator<<(ostream& out, const Graph& g)
{
    if (g._oriented)
        out << "Oriented";
    else out << "Not oriented";
    if (g._weighted)
        out << " weighted";
    out << " graph";
    if (g._flow)
        out <<" with flow";
    out << "\nNumber of vertexes: " << g._nr_vertex
        << "\nNumber of edges: " << g._nr_edge;
    out << "\nEdges are: \n";

    for (int i = 1; i <= g._nr_vertex; ++i)
    {
        for(auto &edge: g._adjacency[i])
        {
             out  << i << " to " << edge.to;
             if (g._weighted)
                out << " with cost " << edge.weight;
             if (g._flow)
                out << " with capacity " << edge.capacity;
            out << "\n";
        }
    }

    return out;
}

Graph::Graph(const bool& oriented, const bool& weighted, const bool& flow)
{
   _nr_vertex = 0;
   _nr_edge = 0;
   _oriented = oriented;
   _weighted = weighted;
   _flow = flow;
}

void Graph::_addEdge(const Edge& e)
{
    /// Add an edge to the graph
    _adjacency[e.from].push_back(Edge(e.from, e.to, e.weight, e.capacity));
    if (!_oriented) _adjacency[e.to].push_back(Edge(e.to, e.from, e.weight, e.capacity));
}

vector<int> Graph::_dist(const int& startV)
{
    /// Solves BFS  returning the distances vector
    //startV = the start vertex of the bfs
    //return: distances = vector of _nr_vertex + 1 length with the distances from the start vertex to the i-th one

    vector<int> distances( _nr_vertex + 1, -1);
    distances[startV] = 0;
    queue<int> order_list;
    order_list.push(startV); // starting BFS from the given vertex
    int current;

    while (!order_list.empty())
    {
        current = order_list.front();
        for (auto &neighbour: _adjacency[current])
        {
            if (distances[neighbour.to] == -1)
            {
                distances[neighbour.to] = distances[current] + 1;
                order_list.push(neighbour.to);
            }
        }
        order_list.pop();
    }
    return distances;
}

int Graph::_bfs(const int& s, int& last)
{
    /// Solves BFS  returning the longest distance
    //s = start node
    // last = last accessed node

    vector<int> distances(_nr_vertex + 1, 0);
    vector<int>visited(_nr_vertex + 1, 0);
    queue<int> order_list;
    order_list.push(s);
    int current, diameter = 0;

    distances[s] = 1;
    visited[s] = 1;
    while (!order_list.empty())
    {
        current = order_list.front();
        order_list.pop();
        for (auto &neighbour: _adjacency[current])
        {
            if (visited[neighbour.to] == 0)
            {
                visited[neighbour.to] = 1;
                distances[neighbour.to] = distances[current] + 1;
                diameter = max(diameter, distances[neighbour.to]); // longest distance starting with s
                order_list.push(neighbour.to);
            }
        }
    }
    last = current;
    return diameter;
}

void Graph::_dfs(const int& start, const int& marker, vector<int>& mark, stack<int>& exit_time)
{
    /// Solves DFS recursively and marks all nodes component
    //start = start vertex
    // marker = value to use when finding a visited vertex
    // mark = vector of _nr_vertex + 1 elements; retains the marker of each vertex
    // exit_time = stack in the order of exiting vertex's DFS

    mark[start] = marker;
    for (auto &neighbour: _adjacency[start])
    {
        if (mark[neighbour.to] == -1)
            _dfs(neighbour.to, marker, mark, exit_time);
    }

    exit_time.push(start);
}

vector<int> Graph::_popEdges(stack<Edge>& st, const Edge& last_edg) // used in leveled_dfs
{
    /// Given an stack of edges pop the elements of the stack till we reach the one given
    // returning the vector of the removed edges
    vector<int> sol;
    int x, y;
    do
    {
        x = st.top().from;
        y = st.top().to;
        st.pop();
        sol.push_back(x);
        sol.push_back(y);
    }while (x != last_edg.from && y !=last_edg.to);

    return sol;
}

void Graph::_popVertex(const int& start, stack<int>& st, int& sol_nr, vector<vector<int>>& sol, vector<bool> &in) //used in leveled_dfs
{
    /// Given a stack of vertex's  pop the elements of the stack till we reach the one given
    // sol_nr = number of strong connected components
    //sol[i] = vertex's of the strong connected component
    //in[i] = the vertex i was visited in the current component

    sol_nr++; // found a new component
    int aux;
    vector<int> aux2;
    do
    {
        aux = st.top();
        in[aux] = 0;
        st.pop();
        aux2.push_back(aux);
    }while (aux != start);
    sol.push_back(aux2); // add the new list to the solution
}

void Graph::_leveledDfs(const int& start, vector<int>& parent, vector<int>& level, vector<int>&  return_level, stack<Edge>& expl_edges, vector<vector<int>>& biconex_comps, vector<vector<int>>& critical_edges)
{
    /// Given a start vertex do a recursive modified DFS, with the purpose of finding the critical vertex's/edges (similar to Tarjan's algorithm)
    // start = current vertex in the dfs
    //level [i] = depth of the vertex in the dfs tree
    //parent[i] = the parent of the vertex i in the dfs tree
    //return_level[i] = the level that the vertex i can return to using return edges
    //expl_edges = stack of edges in the order of discovery - used to discover the biconex components
    //biconex_comps = vector of  the biconex components and their vertex's
    //critical_edges[i] = vector of 2 elements (leetcode restriction) that signify a critical edge

    int nr_children = 0;
    for (auto &child: _adjacency[start])
    {
        nr_children ++; // mark as child
        if (parent[child.to] == -1)
        {
            expl_edges.push(Edge(start, child.to));
            parent[child.to] = start;
            level[child.to] = level[start] + 1;
            return_level[child.to] = level[child.to]; // now the lowest level  reached  from the child is his exact level

            _leveledDfs(child.to, parent, level, return_level, expl_edges, biconex_comps, critical_edges);

            return_level[start] = min(return_level[start], return_level[child.to]); // passed the rest of the DF tree of the child
            // can modify the lowest level reached in case there was a vertex with a return edge

            if (return_level[child.to] > level[start]) // the child cant return higher then his parent
            {
                critical_edges.push_back(vector<int>{start,child.to, 0}); // so this edge is critical in the graph
            }

            if (parent[start] == 0 && nr_children >= 2)// if the root of DF tree has 2 or more  children, means it is a critical vertex
            {
                vector<int> new_comp = _popEdges(expl_edges, {start,child.to});  // get the biconex component of start
                biconex_comps.push_back(new_comp);
            }

            if (parent[start] != 0 && return_level[child.to] >= level[start]) // the child can return to a lower level then his parent
            {
                // so we reached the end of a biconex component
                vector<int> new_comp = _popEdges(expl_edges, {start,child.to}); // get the biconex component of start
                biconex_comps.push_back(new_comp);
            }
        }

        else if (child.to != parent[start]) return_level[start] = min(return_level[start], level[child.to]); // update the lowest level reachable with return edges
    }
}

void Graph::_Tarjan(const int& start, int& time, vector<int>& in_time, vector<int>& return_time, stack<int>& connection, vector<bool>& in_connection, int& nr_ctc, vector<vector<int>>& strong_connected)
{
    /// Given a vertex determine his strong connected component recursively (DF based)
    //start = current vertex
    //time = the nr. of  iterations, meaning the time of discovery of a vertex in DF tree
    //in_time[i] = discovery time of i
    //return_time[i] = the lowest discovery time in the i DF sub-tree  + return edges
    //in_connection[i] = true if i is in the current SCC
    //strong_connected & nr_ctc = nr of strong connected components and the respective vertex's

    time++; // mark a new iteration
    in_time[start] = time;// mark the discovery time of vertex
    return_time[start] = time; // no return edge known
    in_connection[start] = 1;
    connection.push(start);

    for (auto &child: _adjacency[start])
    {
        if (in_time[child.to] == -1)
        {
            _Tarjan(child.to, time,  in_time, return_time, connection, in_connection, nr_ctc, strong_connected);
            return_time[start] = min(return_time[start], return_time[child.to]); // update in case a return edge was found in the child's  df sub-tree
        }
        else if(in_connection[child.to])
        {
            // the child was discovered and in the same SCC as parent => return edge exists
            // update in case child was discovered before the last child that updated
            return_time[start] = min(return_time[start], in_time[child.to]);
       }
    }

    if (return_time[start] == in_time[start]) // a vertex that doesn't have a return edge => end of SCC
    {
        _popVertex(start, connection, nr_ctc, strong_connected, in_connection);
    }
}

bool Graph::_HavelHakimi(const int& n, const int& nr_d, vector<pair<int,int>>& degrees)
{
    /// Given a vector of pairs (Vertex, degree) determine if it is a valid graph
    // and actually makes the graph not just validates
    int  m = nr_d/2; // number of edges is sum of degrees/2
    resize(n);

    vector<vector<bool>> matrix_adj; // access to finding if a edge exists in O(1). but memory + O(n^2)
    matrix_adj.resize(n + 1, vector<bool>(n + 1, 0));

    pair<int,int> max_dg;
    while (m != 0) //the max nr. of iterations is the number of edges in the supposed graph
    {
        sort(degrees.begin(), degrees.end(), SortDescPair); // sort the degrees descending
        if (degrees.size() == 0 || degrees[0].second == 0)
            break; // no vertex left with degree > 0

        max_dg.first = degrees[0].first;
        max_dg.second = degrees[0].second;
        degrees.erase(degrees.begin());

        for (unsigned int k = 0; k < degrees.size() && max_dg.second > 0 ; ++k)
        {
            if (matrix_adj[max_dg.first][degrees[k].first] != 0)
            {
                matrix_adj[max_dg.first][degrees[k].first] = 1; // mark added edge
                matrix_adj[degrees[k].first][max_dg.first] = 1;

                degrees[k].second--; // mark added edge in degrees left
                max_dg.second--;

                if (degrees[k].second < 0) // can't add edges
                    return 0;

                m--;
                AddEdge(max_dg.first, degrees[k].first, 0, 0); // add the edge to the graph
            }
        }

        if (max_dg.second != 0) // didn't add the required edges
            return 0;
    }

    if (m != 0) return 0; // Couldn't add the number of desired edges so not a valid graph

    return 1;
}

int Graph::_flowBfs(const vector<vector<int>>& residual, const int& source, const int& dest, vector<int>& parent, vector<vector<int>>& capacity, vector<vector<int>>& flow)
{
    /// BFS (from start to destination) variation for flux network returning the parent of the destination node given
    queue<int> order_list;
    order_list.push(source);
    parent.assign(_nr_vertex + 1, -1);
    parent[source] = 0;
    int current;

    while (!order_list.empty() && parent[dest] == -1) // we didnt reach the dest
    {
        current = order_list.front();
        order_list.pop();

        // go through the nodes adjacent of the current one
        // and go trough the residual because we update both ways when finding maximum flow
        // maximize the income flux (flow[parent][child])
        // and minimize not used flux (flow[child][parent])
        for (int next: residual[current])
        {
            if (parent[next] == -1 && capacity[current][next] > flow[current][next]) // if it's not visited and can change the flow of the path
            {
                parent[next] = current;
                order_list.push(next);
            }
        }
    }
    return parent[dest]; // return the parent of the dest node (-1 if not found in path)
}

void Graph::_initializeAdMat(vector<vector<int>>& matrix)
{
    // Make copy of adjacency matrix
    matrix.resize(_nr_vertex + 1);

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        for (auto &edge: _adjacency[i])
            matrix[i].push_back(edge.to);
    }
}

int Graph::nrVertex()
{
    return _nr_vertex;
}

int Graph::nrEdges()
{
    return _nr_edge;
}

void Graph::MakeFlowNetwork()
{
    /// Clears adjacency list and makes flow network
    _adjacency.clear();
    _flow = 1;
}

void Graph::MakeNotFlowNetwork()
{
    /// Clears adjacency list and removes flow network marking
    _adjacency.clear();
    _flow = 0;
}

void Graph::MakeOriented()
{
    /// Clears adjacency list and makes oriented graph
    _adjacency.clear();
    _oriented = 1;
}

void Graph::MakeNonOriented()
{
    /// Clears adjacency list and makes non-oriented graph
    _adjacency.clear();
    _oriented = 0;
}

void Graph::MakeWeighted()
{
    /// Clears adjacency list and makes weighted graph
    _adjacency.clear();
    _weighted = 1;
}

void Graph::MakeNonWeighted()
{
    /// Clears adjacency list and makes non-weighted graph
    _adjacency.clear();
    _weighted = 0;
}

void Graph::AddEdge(const int& from, const int& to, const int& cost, const int& capacity)
{
    int realCost = 0, realCap = 0;
    if (_weighted) // Dont allow to add cost if not weighted graph
        realCost = cost;
    if (_flow) // Dont allow to add capacity if not flow network
        realCap = capacity;

    _addEdge(Edge(from, to, realCost, realCap));
    _nr_edge++;
}

void Graph::clear()
{
    _nr_vertex = 0;
    _nr_edge = 0;
    _oriented = 0;
    _weighted = 0;
    _flow = 0;
    _adjacency.clear();
}

void Graph::resize(const int& nr_vertex)
{
     _nr_edge = 0;
    _nr_vertex = nr_vertex;
    _adjacency.clear();
    _adjacency.resize(_nr_vertex + 1);
}

vector<int> Graph::BF(const int& start)
{
    /// Breadth first graph transversal (visiting nodes once)
    // Returns vector with nodes in order
    if(start <1 || start > _nr_vertex)
        return vector<int>();

    vector<int> sol;
    vector<bool> visited(_nr_vertex+ 1, 0);
    queue<int> order_list;
    order_list.push(start);
    sol.push_back(start);
    visited[start] = 1;
    int current;

    while (!order_list.empty())
    {
        current = order_list.front();
        for (auto &neighbour: _adjacency[current])
        {
            if (visited[neighbour.to] == 0)
            {
                sol.push_back(neighbour.to);
                order_list.push(neighbour.to);
                visited[neighbour.to] = 1;
            }
        }
        order_list.pop();
    }

    return sol;
}

vector<int> Graph::DF(const int& start)
{
    /// Depth first graph transversal (visiting nodes once)
    // On the same level the nodes are from level end to level start
    if(start <1 || start > _nr_vertex)
        return vector<int>();

    vector<int> sol;
    stack<int> order;
    vector<int> visited(_nr_vertex + 1, -1);

    _dfs(start, 1, visited, order);

    while (!order.empty())
    {
        sol.push_back(order.top());
        order.pop();
    }

    return sol;
}

pair<int, vector<Edge>> Graph::Prim(const int& start)
{
    /// Prim's Algorithm
    // returning the cost of the MST and the vector of edges that compose it
    if (!_weighted || start < 1 || start > _nr_vertex)
        return make_pair(0, vector<Edge>());

    priority_queue<Edge, vector<Edge>, Edge> heap; // heap of min edge cost
    vector<int> cost(_nr_vertex + 1, INT_MAX); // min cost of every node
    vector<int> parent(_nr_vertex + 1, -1);
    vector<bool> in_APM(_nr_vertex + 1, 0); // keep track of components
    vector<Edge> edges; // for returning the edges in the MST

    int node, total_cost = 0;
    cost[start] = 0; // the cost of staying in place is 0

    heap.push({-1, start, cost[start], 0});

    while (!heap.empty())
    {
        node = heap.top().to;
        heap.pop();
        if (!in_APM[node])
        {
            in_APM[node] = 1;
            total_cost += cost[node]; // it is the min cost  in the heap

            for (auto &neighbour: _adjacency[node]) // go through neighbors
            {
                if (cost[neighbour.to] > neighbour.weight && !in_APM[neighbour.to])
                { // actualize the cost of the current nodes neighbors if their not already in the tree
                    heap.push(neighbour); // add in heap
                    parent[neighbour.to] = node;
                    cost[neighbour.to] = neighbour.weight;
                }
            }

            if (parent[node] != -1) // make sure we dont punt the -1 - start edge that is used to get the cost of  s
                edges.push_back({parent[node], node, cost[node], 0});
        }
    }

    return make_pair(total_cost, edges); // return cost and vector of edges
}

vector<int> Graph::Dijkstra(const int& start)
{
    /// Dijkstra algorithm for finding the minimum cost from start to all the nodes
    // assume there are no negative cycles.
     if (!_weighted || start < 1 || start > _nr_vertex)
        return vector<int>(1,-1);

    priority_queue<Edge, vector<Edge>, Edge> heap;// make heap of min cost
    vector<int> dist(_nr_vertex + 1, INT_MAX);
    vector<bool> in(_nr_vertex + 1, 0);
    int node;
    dist[start] = 0; // cost of staying in place is 0
    in[start] = 1;

    heap.push({0, start, dist[start], 0});
    while (!heap.empty())
    {
        node = heap.top().to;
        heap.pop();

        for (auto &neighbour: _adjacency[node]) //go trough neighbors
        {
             // actualize the cost of getting from start to node
            if (dist[neighbour.to] > dist[node] + neighbour.weight && !in[neighbour.to])
            {
                in[neighbour.to] = 1; // mark visited
                dist[neighbour.to] = dist[node] + neighbour.weight;
                heap.push({neighbour.from, neighbour.to, dist[neighbour.to], 0});
            }
        }
    }

    return dist;
}

int Graph::MaxFlow(const int& source, const int& dest)
{
    /// Edmonds Karp algorithm for computing the maximum flow
    if (!_flow || (source < 1 || source > _nr_vertex) || (dest < 1 || dest > _nr_vertex))
        return -1;

    vector<int> parent;
    vector<vector<int>> residual(_nr_vertex + 1);
    vector<vector<int>> capacity_matrix(_nr_vertex + 1, vector<int>(_nr_vertex + 1, 0)); // for easy access of edge capacity (+ O(n^2) memory and O(m) time for initialization)
    vector<vector<int>> flow_matrix(_nr_vertex + 1, vector<int>(_nr_vertex + 1, 0)); // for easy access of  edge flow

    for (int node = 1; node <= _nr_vertex; ++ node)
    {
        for (auto &edge:_adjacency[node])
        {
            capacity_matrix[node][edge.to] = edge.capacity; // initialize capacity
            residual[node].push_back(edge.to); // mark the neigbours despite the direction (because its residual graph)
            residual[edge.to].push_back(node);
        }
    }

    int max_flow = 0, n_flow;

    while (_flowBfs(residual, source, dest, parent, capacity_matrix, flow_matrix) != -1) // there are optimizations to be made
    {
        for (int k: residual[dest]) // go trough all the nodes adjacent to the destination and optimize those
        {
                if (parent[k] != -1 && capacity_matrix[k][dest] > flow_matrix[k][dest])// only if the edge isnt at max capacity and is found in bfs tree
                {
                       parent[dest] = k; // assume the parent is k
                       n_flow = INT_MAX; // assume the new flow is infinite

                       int  n = dest; // going upwards through the graph's bfs tree
                       while (n != source)
                       {
                           // determine the flow that can be optimized
                           // which is the minimum unused flow of the path
                            n_flow = min(n_flow, capacity_matrix[parent[n]][n] - flow_matrix[parent[n]][n] ); // make the flow of the path
                            n = parent[n];
                       }

                       if (n_flow > 0) // if it actually makes sense to update
                       {
                           n = dest;// going upwards through the graph's bfs tree
                           while (n != source)
                           {
                                flow_matrix[parent[n]][n] += n_flow; //update what the parent sends to child
                                flow_matrix[n][parent[n]] -= n_flow; // update what the child gets from parent
                                n = parent[n];
                           }

                           max_flow += n_flow; // update the maximum flow
                       }
              }
        }
    }

    return max_flow;
}

vector<vector<int>> Graph::RoyFloyd(vector<vector<int>>& dist)
{
    /// Roy-Floyd algorithm for making the minimum cost matrix
    if (!_weighted)
        return (vector<vector<int>>());

    for (int k = 1; k <= _nr_vertex; ++k)
    {
        for (int i = 1; i <= _nr_vertex; ++i)
        {
            for (int j = 1; j <= _nr_vertex; ++j)
            {
                if (dist[i][k] != INT_MAX && dist[k][j] != INT_MAX && dist[i][j] > dist[i][k] + dist[k][j])
                {
                    // update cost from i to j by calculating the costs from an intermediate node
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }

    return dist;
}

int Graph::Diameter()
{
    /// Longest path in a tree
    if (_nr_edge != _nr_vertex - 1)
        return -1;
    // last = last node that gets max path length - so it's the furthest leaf from the root
    int last, aux;
    _bfs(1, last);
    // from last we go back to the root and then descend into the other furthest leaf from root
    return  _bfs(last, aux);
}

vector<int> Graph::DistanceFrom(const int& start)
{
    /// Return distance from start to the other nodes
    if (start < 1 || start > _nr_vertex)
        return vector<int>(1, -1);
    vector<int> result = _dist(start);
    return vector<int> (result.begin() + 1, result.end());
}

int Graph::NrConnectedComponents()
{
    /// Number of connected components of graph
    int result = 0;
    vector<int> components(_nr_vertex + 1, -1);
    stack<int> aux; // not needed in the solving this, but _dfs does multiple things

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (components[i] == -1) // discover new component
        {
            result++;
            _dfs(i, result, components, aux);
        }
    }
    return result;
}

vector<int> Graph::TopologicalSort()
{
    /// Sort DAG topologically
    if (!_oriented)
        return vector<int>();

    vector<int> visited(_nr_vertex + 1, -1);
    vector<int> sol;
    stack<int> aux;

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (visited[i] == -1) // not visited so we go through it
        {
            _dfs(i, 1, visited, aux); // discover all nodes
        }
    }

    while (!aux.empty()) // pop the nodes in the order of exiting DFS
    {
        sol.push_back(aux.top());
        aux.pop();
    }

    return sol; // vector of vertex's in sorted order
}

pair<int,vector<vector<int>>> Graph::BiconnectedComponents()
{
    /// Biconnected components of graph (number and each component)
    vector<vector<int>> sol;
    vector<int> parent(_nr_vertex + 1, -1);
    vector<int> level(_nr_vertex + 1, -1);
    vector<int> rtr_lvl(_nr_vertex + 1, -1);
    stack<Edge> st;
    vector<vector<int>> crt_edg;

    // using a modified DFS
    _leveledDfs(1, parent, level, rtr_lvl, st, sol, crt_edg);

    return pair<int, vector<vector<int>>> (sol.size(), sol); // return (numer of biconex components, sol[i] - vertex's of component)
}

vector<vector<int>> Graph::CriticalConnections()
{
    /// Critical edges of a connected graph
    // With a leveled DF find critical edges
    vector<vector<int>> sol;
    vector<int> parent;
    vector<int> level;
    vector<int> rtr_lvl;
    stack<Edge> st;
    vector<vector<int>> crt_edg;

    _leveledDfs(0, parent, level, rtr_lvl, st, sol, crt_edg);
    return crt_edg; // crt_edg[i] => pair start node - end node of critical edge
}

pair <int, vector<vector<int>>> Graph::StronglyConnectedComponents()
{
    /// Strongly connected components of oriented graph (number and components)
    vector<int> in_time(_nr_vertex + 1, -1);
    vector<int> return_time(_nr_vertex +1, -1);
    stack<int> connection;
    vector<bool> in_connection(_nr_vertex + 1, 0);
    int time;
    time = 0;
    int sol_nr = 0;
    vector<vector<int>> sol;

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (in_time[i] == -1)
            _Tarjan(i, time, in_time, return_time, connection, in_connection, sol_nr, sol);
    }

    return make_pair(sol_nr, sol);
}

vector<int> Graph:: BellmanFord(const int& start)
{
    /// Bellman Ford algorithm for determining the minimum cost from start to the other nodes
    /// and detecting the negative cycles
    // returning vector of costs if there's a negative cycle return is -1
    if (!_weighted || start < 1 || start > _nr_vertex)
       return vector<int>(1,-1);

    vector<int> dist(_nr_vertex + 1, INT_MAX); // assume we cant get to any node
    vector<int> parent(_nr_vertex + 1, -1);
    queue<int> order;
    vector<bool> in_q(_nr_vertex + 1, 0);
    vector<int> count_in_q(_nr_vertex  + 1, 0); // how many times was N in queue

    int node;
    dist[start] = 0;
    in_q[start] = 1; // mark that the current node is in queue
    count_in_q[start] = 1;
    order.push(start);

    while(!order.empty())
    {
       node = order.front();
       order.pop();
       in_q[node] = 0; // mark that it was popped out of the queue

        for (auto &neighbour:_adjacency[node]) // go through neighbors
        {
            if (dist[neighbour.to] > dist[node] + neighbour.weight) // we can actualize the distance from node to neighbor
            {
                dist[neighbour.to] = dist[node] + neighbour.weight;
                parent[neighbour.to] = node; // mark from who we got to neighbor
                if (!in_q[neighbour.to])
                {
                     // if the node was in queue more times then nodes there is a negative cycle
                     if (count_in_q[neighbour.to] > _nr_vertex)
                        return vector<int>(1,-1);
                    order.push(neighbour.to);
                    count_in_q[neighbour.to] += 1;
                }
            }
        }
    }

    return dist;
}

tuple<int, int, vector<Edge>> Graph::MinimumSpanningTree()
{
    ///Using Prim's Algorithm return the cost of MST, number of edges (nr vertex - 1) and the actual edges
    if (!_weighted)
        return make_tuple(-1,-1,vector<Edge>());
    pair<int, vector<Edge>> solution = Prim(1);
    return make_tuple(solution.first, solution.second.size(), solution.second);
}

void Graph::HavelHakimi(istream &in)
{
    /// Solving the problem "determine if array of degrees is a valid non-oriented graph" with Havel Hakimi theorem
    /// and making the graph. Destroys current graph.
    if (_oriented)
        return;

    bool breakflag = 0; // flag in case the algorithm should end so we can display first
    int aux, n, m = 0, nr_zeros = 0;
    vector<pair<int,int>> degrees; // <node, degree>

    in >> n;
    for (int i = 1; i <= n; ++i)
    {
        in >> aux;
        // if degree > number of vertex's => impossible
        if (aux > n)
        {
            breakflag = 1;
            break;
        }
        if (aux == 0)
            nr_zeros++;

        m += aux;
        degrees.push_back(make_pair(i, aux));
    }
    if  (m % 2 == 1 || breakflag )
        cout << "Not a graph.\n";
    else if (nr_zeros == n)
    {
        resize(n);
        cout <<"An empty graph with " << n << " vertex's.\n";
    }
    else
     {
            bool answer = _HavelHakimi(n,m,degrees);
            if  (!answer)
            {
                clear();
                cout << "Not a graph.\n";
            }
            else
            {
                cout << *this;
            }
     }
}

vector<int> Graph::EulerianCycle()
{
    /// Determine Eulerian cycle of graph and return the nodes component
    /// if there isn't any return -1
    // Reunion of disjunct cycles

    //make another adjacency list not to destroy current one (+O(m) time and O(n^2) memory)
    // but it is more reasonable that the graph is not destroyed when calling a method.
    vector<vector<int>> adj_list;
    _initializeAdMat(adj_list);

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (adj_list[i].size() % 2 == 1) // theorem
            return vector<int>(1, -1);
    }

    stack<int> nodes; // stack of nodes to make recursive alg. iterative
    nodes.push(1); // start from whatever node (cycle theory)

    vector<int> cycle; // solution
    int current, next_node, nr = 0;

    while (!nodes.empty())
    {
        current = nodes.top(); // take the current node

        if (nr == _nr_edge) // got Eulerian cycle so break loop
            break;
        if (adj_list[current].empty()) // passed through all neighbors
        {
            nr++; // mark nr of edges used
            nodes.pop();
            cycle.push_back(current); // add node to cycle
        }
        else
        {
                next_node = adj_list[current].back(); // get a random neighbor
                adj_list[current].pop_back(); // remove the edge from current node

                auto pos = find(adj_list[next_node].begin(), adj_list[next_node].end(), current);
                if (pos != adj_list[next_node].end())
                    adj_list[next_node].erase(pos); // remove the edge from the other node side
                nodes.push(next_node); // add node to recursion
        }
    }

    if (nr != _nr_edge) // there's no cycle.
        return vector<int> (1, -1);

    return cycle;
}

bool Graph::VerifyHavelHakimi(const int& n, const int& nr_d, vector<pair<int,int>>& degrees)
{
    /// Just verifies the theorem
    int  m = nr_d/2; // number of edges is sum of degrees/2

    vector<vector<bool>> matrix_adj; // access to finding if a edge exists in O(1). but memory + O(n^2)
    matrix_adj.resize(n + 1, vector<bool>(n + 1, 0));

    pair<int,int> max_dg;
    while (m != 0) //the max nr. of iterations is the number of edges in the supposed graph
    {
        sort(degrees.begin(), degrees.end(), SortDescPair); // sort the degrees descending
        if (degrees.size() == 0 || degrees[0].second == 0)
            break; // no vertex left with degree > 0

        max_dg.first = degrees[0].first;
        max_dg.second = degrees[0].second;
        degrees.erase(degrees.begin());

        for (unsigned int k = 0; k < degrees.size() && max_dg.second > 0 ; ++k)
        {
            if (matrix_adj[max_dg.first][degrees[k].first] != 0)
            {
                matrix_adj[max_dg.first][degrees[k].first] = 1; // mark added edge
                matrix_adj[degrees[k].first][max_dg.first] = 1;

                degrees[k].second--; // mark added edge in degrees left
                max_dg.second--;

                if (degrees[k].second < 0) // can't add edges
                    return 0;

                m--;
            }
        }

        if (max_dg.second != 0) // didn't add the required edges
            return 0;
    }

    if (m != 0) return 0; // Couldn't add the number of desired edges so not a valid graph

    return 1;
}

/// Pentru meniu
class Infoarena
{
protected:
    Infoarena()
    {
    }

    static Infoarena* instance;

    bool check(ifstream& in)
    {
        if (!in.good())
        {
            cout << "Verificati daca fisierul de intrare exista.\n";
            in.close();
            return false;
        }
        return true;
    }

    void showRead(const char* filename)
    {
        cout << "Citire date din \""<< filename <<".in\".\n";
    }

    void showEnd(const char* filename)
    {
        cout << "Verificati fisierul de iesire \"" << filename <<".out\".\n";
    }

public:
    Infoarena(Infoarena&) = delete;
    void operator=(const Infoarena&) = delete;

    static Infoarena* init()
    {
        if (instance == NULL)
            instance = new Infoarena();
        return instance;
    }

    void havelhakimi()
    {
        Graph g;
        g.HavelHakimi(cin);
    }

    void bfs()
    {
        showRead("bfs");
        ifstream in("bfs.in");
        if (!check(in))
            return;

        int n, m, s, x, y;
        in >> n >> m >> s;
        Graph g(1);
        g.resize(n);
        for (int i = 1; i <=m; ++i)
        {
            in >> x >> y;
            g.AddEdge(x,y);
        }
        vector<int> sol = g.DistanceFrom(s);
        ofstream out("bfs.out");
        for (unsigned int i = 0; i < sol.size(); ++i)
            out << sol[i] << " ";
        out.close();
        showEnd("bfs");
    }

    void dfs()
    {
        showRead("dfs");
        ifstream in("dfs.in");
        if (!check(in))
            return;
        Graph g;
        in >> g;
        in.close();
        int sol = g.NrConnectedComponents();

        ofstream out("dfs.out");
        out << sol;
        out.close();
        showEnd("dfs");
    }

    void sortaret()
    {
        showRead("sortaret");
        ifstream in("sortaret.in");
                if (!check(in))
            return;

        Graph g(1);
        in >> g;
        in.close();

        vector<int> sol = g.TopologicalSort();
        ofstream out("sortaret.out");
        for (unsigned int i = 0; i < sol.size(); ++i)
            out << sol[i] << " ";
        out.close();
        showEnd("sortaret");
    }

    void biconex()
    {
        showRead("biconex");
        ifstream in("biconex.in");
        if (!check(in))
            return;

        Graph g;
        in >> g;
        in.close();

        ofstream out("biconex.out");
        pair<int, vector<vector<int>>> sol = g.BiconnectedComponents();
        out << sol.first << "\n";

        for (int i = 0; i < sol.first; ++i)
        {
            sort(sol.second[i].begin(), sol.second[i].end());
            sol.second[i].erase(unique(sol.second[i].begin(), sol.second[i].end()), sol.second[i].end());
            for (unsigned int e_idx = 0; e_idx < sol.second[i].size(); ++e_idx)
            {
                out << sol.second[i][e_idx] << " ";
            }
            out << "\n";
        }
        out.close();
        showEnd("biconex");
    }

    void ctc()
    {
        showRead("ctc");
        ifstream in("ctc.in");
        if (!check(in))
            return;
        Graph g(1);
        in >> g;
        in.close();
        pair<int, vector<vector<int>>> solution = g.StronglyConnectedComponents();
        ofstream out("ctc.out");
        out << solution.first << endl;
        for (int i = 0; i < solution.first; ++i)
        {
            for (auto &it: solution.second[i])
                out << it << " ";
            out << "\n";
        }
        out.close();
        showEnd("ctc");
    }

    void leetcc()
    {
        showRead("criticalconnections");
        ifstream in("criticalconnections.in");
         if (!check(in))
            return;
        ofstream out("criticalconnections.out");
        Graph g;
        in >> g;
        in.close();
        vector<vector<int>> crt_edg = g.CriticalConnections();

        for (unsigned int i = 0; i < crt_edg.size(); ++i)
        {
            for (unsigned int e_idx = 0; e_idx < crt_edg[i].size(); ++e_idx)
            {
                out << crt_edg[i][e_idx] << " ";
            }
            out << "\n";
        }
        out.close();
        showEnd("criticalconnections");
    }

    void apm()
    {
        showRead("apm");
        ifstream in("apm.in");
         if (!check(in))
            return;
        Graph g(0,1);
        in >> g;
        in.close();
        ofstream out("apm.out");
        tuple<int,int, vector<Edge>> sol = g.MinimumSpanningTree();
        out << get<0>(sol) << "\n" << get<1>(sol) << "\n";

        for (auto &edge:get<2>(sol))
        {
            out << edge.from << " " << edge.to << "\n";
        }

        out.close();
        showEnd("apm");
    }

    void dijkstra()
    {
        showRead("dijkstra");
        ifstream in("dijkstra.in");
         if (!check(in))
            return;
        Graph g(1,1);
        in >> g;
        in.close();
        vector<int> solution = g.Dijkstra(1);
        ofstream out("dijkstra.out");
        int n = g.nrVertex();
        for (int i = 2;  i < n + 1 ; ++i)
        {
            if (solution[i] == INT_MAX) out << 0 << " ";
            else out << solution[i] << " ";
        }
        out.close();
         showEnd("dijkstra");
    }

    void disjoint()
    {
        showRead("disjoint");
        ifstream in("disjoint.in");
         if (!check(in))
            return;
        int n, m;
        in >> n >> m;

        vector<Forest> f = Forest::MakeForest(n);
        int x, y, cod;

        ofstream out("disjoint.out");
        while (m)
        {
            m--;

            in >> cod >> x >> y;
            if (cod == 1)   Forest::UnionForest(x, y, f);
            if (cod == 2)
            {
                if (Forest::CheckForest(x,y,f) == 1)
                    out << "DA\n";
                else out  << "NU\n";
            }
        }
        out.close();
        showEnd("disjoint");
    }

    void bellmanford()
    {
        showRead("bellmanford");
        ifstream in("bellmanford.in");
        if (!check(in))
            return;
        Graph g(1,1);
        in >> g;
        in.close();
        ofstream out("bellmanford.out");
        vector<int> sol = g.BellmanFord(1);
        if (sol.size() == 1 && sol[0] == -1)
            out << "Ciclu negativ!\n";
        else
        {
            for (int i = 2; i <= g.nrVertex(); ++i)
                out << sol[i] << " ";
        }
        out.close();
        showEnd("bellmanford");
    }

    void maxflow()
    {
        showRead("maxflow");
        ifstream in("maxflow.in");
        if (!check(in))
            return;
        Graph g(1, 0, 1);
        in >> g;
        in.close();
        ofstream out("maxflow.out");
        out << g.MaxFlow(1, g.nrVertex());
        out.close();
        showEnd("maxflow");
    }

    void royfloyd()
    {
        showRead("royfloyd");
        ifstream in("royfloyd.in");
        if (!check(in))
            return;
        int n;
        in >> n;
        Graph g(1, 1);
        g.resize(n);
        vector<vector<int>> sol (n + 1, vector<int> (n+1, INT_MAX));
        for (int i = 1; i <= n; ++i)
        {
            for (int j = 1; j <= n; ++j)
            {
                in >> sol[i][j];
                if (sol[i][j] == 0 && i != j)
                    sol[i][j] = INT_MAX;
            }
        }
        in.close();
        g.RoyFloyd(sol);
        ofstream out("royfloyd.out");
        for (int i = 1; i <= n; ++i)
        {
            for (int j = 1; j <= n; ++j)
            {
                if (sol[i][j] == INT_MAX)
                    out << 0 << " ";
                else out << sol[i][j] << " ";
            }
            out <<"\n";
        }
        out.close();
        showEnd("royfloyd");
    }

    void darb()
    {
        showRead("darb");
        ifstream in("darb.in");
        if (!check(in))
            return;
        int n, x, y;
        Graph g;
        in >> n;
        g.resize(n);
        for (int i = 1; i <= n - 1; ++i)
        {
            in >> x >> y;
            g.AddEdge(x,y);
        }
        in.close();
        ofstream out("darb.out");
        out << g.Diameter();
        out.close();
        showEnd("darb");
    }

    void ciclueuler()
    {
        showRead("ciclueuler");
        ifstream in("ciclueuler.in");
        if (!check(in))
            return;
        Graph g;
        in >> g;
        in.close();
        ofstream out("ciclueuler.out");
        vector<int> sol = g.EulerianCycle();
        if (sol.size() == 1 && sol[0] == -1)
            out << "-1\n";
        else
        {
            for (auto &node: sol)
                out << node << " ";
        }
        out.close();
        showEnd("ciclueuler");
    }
};

Infoarena *Infoarena::instance = NULL;
Infoarena *probleme = probleme->init();

void displayOption(const int& i, const char* option)
{
    cout << i << ")\tProblema " << option << "\n";
}

void meniu()
{
     int code = 0;

    for (int i = 1; i <= 16; ++i)
    {
        switch (i)
        {
            case 16: cout <<"0 for EXIT\n"; break;
            case 1: cout << "Tema \"BF-DF si aplicatii\":\n"; displayOption(i, "bfs"); break;
            case 2:  displayOption(i, "dfs"); break;
            case 3:  displayOption(i, "sortare topologica"); break;
            case 4:  displayOption(i, "biconex"); break;
            case 5:  displayOption(i, "ctc"); break;
            case 6:  displayOption(i, "Havel Hakimi"); break;
            case 7: displayOption(i, "critical connections (leetcode)"); break;
            case 8: cout << "\nTema \"Drumuri minime si APM\":\n"; displayOption(i, "apm"); break;
            case 9: displayOption(i, "Dijkstra"); break;
            case 10: displayOption(i, "disjoint"); break;
            case 11: displayOption(i, "Bellman-Ford"); break;
            case 12: cout << "\nTema pentru laboratorul 5\n"; displayOption(i, "maxflow"); break;
            case 13:  displayOption(i, "Roy-Floyd"); break;
            case 14: displayOption(i, "darb"); break;
            case 15: cout << "\nTema pentru laboratorul 6\n"; displayOption(i, "Ciclu Eulerian"); break;
        };
    }

    do
    {
        cout << "Optiune: ";
        cin >> code;
        switch (code)
        {
            case 1: probleme->bfs(); break;
            case 2: probleme->dfs(); break;
            case 3: probleme->sortaret(); break;
            case 4: probleme->biconex(); break;
            case 5: probleme->ctc(); break;
            case 6: probleme->havelhakimi(); break;
            case 7: probleme->leetcc(); break;
            case 8: probleme->apm(); break;
            case 9: probleme->dijkstra(); break;
            case 10: probleme->disjoint(); break;
            case 11: probleme->bellmanford(); break;
            case 12: probleme->maxflow(); break;
            case 13: probleme->royfloyd(); break;
            case 14: probleme->darb(); break;
            case 15: probleme->ciclueuler(); break;
            case 0: return; break;
            default: cout << "Optiune invalida.\n";
        };
    } while(code != 0);
}

int main()
{
    // Se poate accesa meniu pentru lista de probleme rezolvate
    // la teme. Sau a se scrie un program propriu ce foloseste
    // clasa Graph.

    Graph g(1);
    cin >> g;
    vector<int> bfs = g.BF(1);

    for (auto elem: bfs)
        cout << elem << " ";
    cout << endl;

    vector<int> dfs = g.DF(1);

    for (auto elem: dfs)
        cout << elem << " ";
    cout << endl;

    cout << g.MaxFlow(1,4) << endl;

    g.MakeFlowNetwork();
    cin >> g;

    cout << g.MaxFlow(1, 4) << endl;

    return 0;
}
