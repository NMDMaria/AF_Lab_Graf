#include <bits/stdc++.h>

using namespace std;

class Edge
{
    /*
    class that helps with the definition of edges, with costs and or capacity
    */
private:
    int from, where;
    int weight;
    int capacity;
public:
    friend class Graph;
    friend void infoarenaAPM();
    Edge(const int& f=0, const int& w=0, const int& wg = 0, const int& cap = 0)
    {
        from = f;
        where = w;
        weight = wg;
        capacity = cap;
    }
     bool operator()(const Edge &left, const Edge& right)
    {
        return left.weight > right.weight;
    }
};


bool pairsortsnddesc(const pair<int,int>& i, const pair<int,int>& j)
{
     /// Function used to sort a vector of pairs descending by the second component
    return i.second > j.second;
}

struct Forest
{
    int parent, size;
};


class Graph
{
 private:
    int _nr_vertex;
    int _nr_edge;
    bool _oriented;
    bool _weighted;
    bool _flow;
    vector<vector<Edge>> _adjacency;

    // Internal
    vector<int> _popEdges(stack<Edge>& , const Edge&);
    void _popVertex(const int&, stack<int>&, int&, vector<vector<int>>&, vector<bool>&);
    int  _bfs(const int&, int&);
    vector<int> _bfs(const int&);
    void _dfs(const int&, const int&, vector<int>&, stack<int>&);
    void _leveled_dfs(const int&, vector<int>& , vector<int>&, vector<int>& , stack<Edge>&, vector<vector<int>>&, vector<vector<int>>&);
    void _tarjan(const int&, int&, vector<int>&, vector<int>&, stack<int>&, vector<bool>& , int&, vector<vector<int>>&);
    bool _HavelHakimi(const int&, const int&, vector<pair<int,int>>& );
    void _addEdge(const Edge&);
    void _resize(const int&, const int&, const bool&, const bool&, const bool&);
    vector<int> _BellmanFord(const int&);
    vector<Forest> _makeForest(const int&);
    int _findParent(int, vector<Forest>&);
    void _unionForest(int , int , vector<Forest>& );
    bool _checkForest(int, int, vector<Forest>&);
    int _flowBfs(const vector<vector<int>>&, const int& source, const int& dest, vector<int>& parent, vector<vector<int>>& capacity, vector<vector<int>>& flow);
    void _initializeAdMat(vector<vector<int>>&);

 public:
    Graph(const int& n = 0, const int& m = 0, const bool& o = 0, const bool& w = 0, const bool &flow = 0);
    ifstream& readEdges(ifstream&);

    pair<int,vector<Edge>> Prim();
    vector<vector<int>> RoyFloyd(ifstream& in);
    vector<int> BFS(const int&);
    int CompConexe();
    vector<int> sortareTopologica();
    pair<int,vector<vector<int>>> compBiconexe();
    vector<vector<int>> criticalConnections();
    pair<int, vector<vector<int>>> componenteTareConexe();
    void HavelHakimi(ifstream&);
    tuple<int, int, vector<Edge>> arborePartialCostMinim();
    vector<int> Dijkstra(const int&);
    vector<int> BellmanFord();
    int  maxFlow(const int&, const int &);
    vector<bool> paduriDisjuncte(int, int, ifstream &);
    int Edmonds_Karp();
    int Darb(ifstream& in);
    vector<int> EulerianCycle();
};

Graph::Graph(const int& n, const int& m, const bool& o, const bool& w, const bool& f)
{
   _nr_vertex= n;
   _nr_edge = m;
   _oriented = o;
   _weighted = w;
   _flow = f;
   if (n != 0) _adjacency.resize(n + 1);
}

ifstream& Graph::readEdges(ifstream& in)
{
    /// Reads the edges of the graph
    int x, y, c = 0, capacity = 0;
    for(int i = 0; i < _nr_edge; ++i)
    {
        in >> x >> y;
        if (_weighted == 1)
        {
            in >> c;
        }
        if (_flow == 1)
        {
            in >> capacity;
        }
        _addEdge({x, y, c, capacity});
    }
    return in;
}

/// Definitions of internal methods
void Graph::_resize(const int& new_n, const int& new_m, const bool& new_o, const bool& new_w, const bool& new_f)
{

    _nr_vertex= new_n;
    _nr_edge = new_m;
    _oriented = new_o;
    _weighted = new_w;
    _flow = new_f;
    _adjacency.clear();
    _adjacency.resize(new_n + 1);
}

void Graph::_addEdge(const Edge& e)
{
    /// Add an edge to the graph
    _adjacency[e.from].push_back(Edge(e.from, e.where, e.weight, e.capacity));
    if (!_oriented) _adjacency[e.where].push_back(Edge(e.where, e.from, e.weight, e.capacity));
}

vector<Forest> Graph::_makeForest(const int& N)
{
    vector<Forest> F(N + 1);
    for (int i = 1; i <= N; ++i)
    {
        F[i].parent = i;
        F[i].size = 1;
    }

    return F;
}

int Graph::_findParent(int x, vector<Forest>& F)
{
    while (F[x].parent != x)
    {
        F[x].parent = F[ F[x].parent].parent;
        F[x].size = F[F[x].parent].size;
        x = F[x].parent;
    }
    return x;
}

void Graph::_unionForest(int x, int y, vector<Forest>& F)
{
    x = _findParent(x, F);
    y = _findParent(y, F);

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

bool Graph::_checkForest(int x, int y, vector<Forest>& F)
{
    if (_findParent(x, F) == _findParent(y, F))
        return 1;
    return 0;
}

vector<int> Graph::_bfs(const int& startV)
{
    /// Solves BFS  returning the distances vector
    /*startV = the start vertex of the bfs
     return: distances = vector of _nr_vertex + 1 length with the distances from the start vertex to the i-th one
    */
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
            if (distances[neighbour.where] == -1)
            {
                distances[neighbour.where] = distances[current] + 1;
                order_list.push(neighbour.where);
            }
        }
        order_list.pop();
    }
    return distances;
}

int Graph::_bfs(const int& s, int& last)
{
    /// Solves BFS  returning the longest distance
    /*startV = the start vertex of the bfs
     return: distances = vector of _nr_vertex + 1 length with the distances from the start vertex to the i-th one
    */

    vector<int> distances(_nr_vertex + 1, 0);
    vector<int>visited(_nr_vertex + 1, 0);
    queue<int> order_list;
    order_list.push(s);
    int current, diametru = 0;

    distances[s] = 1;
    visited[s] = 1;
    while (!order_list.empty())
    {
        current = order_list.front();
        order_list.pop();
        for (auto &neighbour: _adjacency[current])
        {
            if (visited[neighbour.where] == 0)
            {
                visited[neighbour.where] = 1;
                distances[neighbour.where] = distances[current] + 1;
                diametru = max(diametru, distances[neighbour.where]);
                order_list.push(neighbour.where);
            }
        }
    }
    //cout << current << " " <<  distances[current] << endl;
    last = current;
    return diametru;
}

void Graph::_dfs(const int& start, const int& marker, vector<int>& mark, stack<int>& exit_time)
{
    /// Solves DFS recursively
    /*
     start = start vertex
     marker = value to use when finding a visited vertex
     mark = vector of _nr_vertex + 1 elements; retains the marker of each vertex
     exit_time = stack in the order of exiting vertex's DFS
    */

    mark[start] = marker;
    for (auto &neighbour: _adjacency[start])
    {
        if (mark[neighbour.where] == -1)
            _dfs(neighbour.where, marker, mark, exit_time);
    }

    exit_time.push(start);
}

vector<int> Graph::_popEdges(stack<Edge>& st, const Edge& last_edg)
{
    /// Given an stack of edges pop the elements of the stack till we reach the one given
    // returning the vector of the removed edges
    vector<int> sol;
    int x, y;
    do
    {
        x = st.top().from;
        y = st.top().where;
        st.pop();
        sol.push_back(x);
        sol.push_back(y);
    }while (x != last_edg.from && y !=last_edg.where);

    return sol;
}

void Graph::_popVertex(const int& start, stack<int>& st, int& sol_nr, vector<vector<int>>& sol, vector<bool> &in)
{
    /// Given a stack of vertex's  pop the elements of the stack till we reach the one given
    /* sol_nr = number of strong connected components
     sol[i] = vertex's of the strong connected component
     in[i] = the vertex i was visited in the current component
    */
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

void Graph::_leveled_dfs(const int& start, vector<int>& parent, vector<int>& level, vector<int>&  return_level, stack<Edge>& expl_edges, vector<vector<int>>& biconex_comps, vector<vector<int>>& critical_edges)
{
    /// Given a start vertex do a recursive modified DFS, with the purpose of finding the critical vertex's/edges (similar to Tarjan's algorithm)
    /* start = current vertex in the dfs
     level [i] = depth of the vertex in the dfs tree
     parent[i] = the parent of the vertex i in the dfs tree
     return_level[i] = the level that the vertex i can return to using return edges
     expl_edges = stack of edges in the order of discovery - used to discover the biconex components
     biconex_comps = vector of  the biconex components and their vertex's
      critical_edges[i] = vector of 2 elements (leetcode restriction) that signify a critical edge
    */
    int nr_children = 0;
    for (auto &child: _adjacency[start])
    {
        nr_children ++; // mark as child
        if (parent[child.where] == -1)
        {
            expl_edges.push(Edge(start, child.where));
            parent[child.where] = start;
            level[child.where] = level[start] + 1;
            return_level[child.where] = level[child.where]; // now the lowest level  reached  from the child is his exact level

            _leveled_dfs(child.where, parent, level, return_level, expl_edges, biconex_comps, critical_edges);

            return_level[start] = min(return_level[start], return_level[child.where]); // passed the rest of the DF tree of the child
            // can modify the lowest level reached in case there was a vertex with a return edge

            if (return_level[child.where] > level[start]) // the child cant return higher then his parent
            {
                critical_edges.push_back(vector<int>{start,child.where, 0}); // so this edge is critical in the graph
            }

            if (parent[start] == 0 && nr_children >= 2)// if the root of DF tree has 2 or more  children, means it is a critical vertex
            {
                vector<int> new_comp = _popEdges(expl_edges, {start,child.where});  // get the biconex component of start
                biconex_comps.push_back(new_comp);
            }

            if (parent[start] != 0 && return_level[child.where] >= level[start]) // the child can return to a lower level then his parent
            {
                // so we reached the end of a biconex component
                vector<int> new_comp = _popEdges(expl_edges, {start,child.where}); // get the biconex component of start
                biconex_comps.push_back(new_comp);
            }
        }

        else if (child.where != parent[start]) return_level[start] = min(return_level[start], level[child.where]); // update the lowest level reachable with return edges
    }
}

void Graph::_tarjan(const int& start, int& time, vector<int>& in_time, vector<int>& time_return_vert, stack<int>& connection, vector<bool>& in_connection, int& nr_ctc, vector<vector<int>>& strong_connected)
{
    /// Given a vertex determine his strong connected component recursively (DF based)
    /* start = current vertex
     time = the nr. of  iterations, meaning the time of discovery of a vertex in DF tree
     in_time[i] = discovery time of i
     time_return_vert[i] = the lowest discovery time in the i DF sub-tree  + return edges
     in_connection[i] = true if i is in the current SCC
     strong_connected & nr_ctc = nr of strong connected components and the respective vertex's
    */
    time++; // mark a new iteration
    in_time[start] = time;// mark the discovery time of vertex
    time_return_vert[start] = time; // no return edge known
    in_connection[start] = 1;
    connection.push(start);

    for (auto &child: _adjacency[start])
    {
        if (in_time[child.where] == -1)
        {            // continue df
            _tarjan(child.where, time,  in_time, time_return_vert, connection, in_connection, nr_ctc, strong_connected);
            time_return_vert[start] = min(time_return_vert[start], time_return_vert[child.where]); // update in case a return edge was found in the child's  df sub-tree
        }
        else if(in_connection[child.where]) time_return_vert[start] = min(time_return_vert[start], in_time[child.where]);
        // the child was discovered and in the same SCC as parent => return edge exists
       // update in case child was discovered before the last child that updated
    }

    if (time_return_vert[start] == in_time[start]) // a vertex that doesn't have a return edge => end of SCC
    {
        _popVertex(start, connection, nr_ctc, strong_connected, in_connection);
    }
}

bool Graph::_HavelHakimi(const int& n, const int& nr_d, vector<pair<int,int>>& degrees)
{
    /// Given a vector of pairs (Vertex, degree) determine if it is a valid graph
    int  m = nr_d/2; // number of edges is sum of degrees/2
    _resize(n, m, 0, 0, 0);

    vector<vector<bool>> matrix_adj; // access to finding if a edge exists in O(1). but memory + O(n^2)
    matrix_adj.resize(n + 1);
    for (int i = 1; i <= n; ++i)
        matrix_adj[i].resize(n + 1, 0);
    while (m) //the max nr. of iterations is the number of edges in the supposed graph
    {
        sort(degrees.begin(), degrees.end(), pairsortsnddesc); // Sort the degrees descending
        if (degrees.size() == 0 || degrees[0].second == 0) break; // No vertex left with degree > 0
        pair<int,int> max_dg;
        max_dg.first = degrees[0].first;
        max_dg.second = degrees[0].second;
        degrees.erase(degrees.begin());
        for (unsigned int k = 0; k < degrees.size() && max_dg.second > 0 ; ++k)
        {
            if (!matrix_adj[max_dg.first][degrees[k].first])
            {
                matrix_adj[max_dg.first][degrees[k].first] = matrix_adj[degrees[k].first][max_dg.first] = 1;
                degrees[k].second--;
                max_dg.second--;
                if (degrees[k].second < 0) return 0;
                m--;

                _addEdge({max_dg.first, degrees[k].first, 0, 0});
            }
        }
        if (max_dg.second != 0) return 0;
    }

    if (m!=0) return 0; // Couldn't add the number of desired edges so not a valid graph
    return 1;
}

vector<int> Graph:: _BellmanFord(const int& start)
{
    vector<int> dist(_nr_vertex + 1, INT_MAX); // assume we cant get to any node
    vector<int> parent(_nr_vertex + 1, -1);
    queue<int> order;
    vector<bool> in_q(_nr_vertex + 1, 0);// mark if we have
    vector<int> count_in_q(_nr_vertex  + 1, 0); // how many times was N in queue

    int node;
    dist[start] = 0;
    in_q[start] = 1;
    count_in_q[start] = 1;
    order.push(start);
    // based on BF
    while(!order.empty())
    {
       node = order.front();
       order.pop();
       in_q[node] = 0; // mark we don't have it rn in the queue

        for (auto &neighbour:_adjacency[node])
        {
            if (dist[neighbour.where] > dist[node] + neighbour.weight)
            { // we can actualize the distance from node to neighbor
                dist[neighbour.where] = dist[node] + neighbour.weight;
                parent[neighbour.where] = node; // mark from who we got to neighbor
                if (!in_q[neighbour.where])
                {
                     if (count_in_q[neighbour.where] > _nr_vertex) return vector<int>(1,-1);
                     // if the node was in queue more times then nodes there is a negative cycle
                    order.push(neighbour.where);
                    count_in_q[neighbour.where] += 1;
                }
            }
        }
    }

    return dist;
}

int Graph::_flowBfs(const vector<vector<int>>& residual, const int& source, const int& dest, vector<int>& parent, vector<vector<int>>& capacity, vector<vector<int>>& flow)
{
    // Bfs for capacity graph from source to dest
    queue<int> order_list;
    order_list.push(source);
    parent.assign(_nr_vertex + 1, -1);
    parent[source] = 0;
    int current;

    while (!order_list.empty() && parent[dest] == -1) // we didnt reach the dest
    {
        current = order_list.front();
        order_list.pop();

        for (int next: residual[current]) // we go through all the nodes adjacent to current node
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


///Definitions of public functions
pair<int, vector<Edge>> Graph::Prim()
{
    /// Prim's Algorithm
    priority_queue<Edge, vector<Edge>, Edge> heap; // heap of min edge cost
    vector<int> cost(_nr_vertex + 1, INT_MAX); // min cost of every node
    vector<int> parent(_nr_vertex + 1, -1);
    vector<bool> in_APM(_nr_vertex + 1, 0); // keep track of components
    vector<Edge> edges; // for returning the edges in the MST

    int node, total_cost = 0;
    cost[1] = 0; // the cost of staying in place is 0

    heap.push({-1, 1, cost[1], 0});

    while (!heap.empty())
    {
        node = heap.top().where;
        heap.pop();
        if (!in_APM[node])
        {
            in_APM[node] = 1;
            total_cost += cost[node]; // it is the min cost  in the heap

            for (auto &neighbour: _adjacency[node]) // go through neighbors
            {
                if (cost[neighbour.where] > neighbour.weight && !in_APM[neighbour.where])
                { // actualize the cost of the current nodes neighbors if their not already in the tree
                    heap.push(neighbour); // add in heap
                    parent[neighbour.where] = node;
                    cost[neighbour.where] = neighbour.weight;
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
    priority_queue<Edge, vector<Edge>, Edge> heap;// make heap of min cost
    vector<int> dist(_nr_vertex + 1, INT_MAX);
    vector<bool> in(_nr_vertex + 1, 0);
    int node;
    dist[start] = 0; // cost of staying in place is 0

    heap.push({0, start, dist[start], 0});
    while (!heap.empty())
    {
        node = heap.top().where;
        heap.pop();

        for (auto &neighbour: _adjacency[node]) //go trough neighbors
        {
            if (dist[neighbour.where] > dist[node] + neighbour.weight && !in[neighbour.where])
            {// we can actualize the cost of getting from start to node
                in[neighbour.where] = 1; // mark visited
                dist[neighbour.where] = dist[node] + neighbour.weight;
                heap.push({neighbour.from, neighbour.where, dist[neighbour.where], 0});
            }
        }
    }

    return vector<int>(dist.begin() + 2, dist.end()); // for this we don't need the dist[start]
}

int Graph::maxFlow(const int& source, const int& dest)
{
    vector<int> parent;
    vector<vector<int>> residual(_nr_vertex + 1);
    vector<vector<int>> capacity_matrix(_nr_vertex + 1, vector<int>(_nr_vertex + 1, 0)); // for easy access of edge capacity
    vector<vector<int>> flow_matrix(_nr_vertex + 1, vector<int>(_nr_vertex + 1, 0)); // for easy access of  edge flow

    for (int node = 1; node <= _nr_vertex; ++ node)
    {
        for (auto &edge:_adjacency[node])
        {
            capacity_matrix[node][edge.where] = edge.capacity; // initialize capacity
            residual[node].push_back(edge.where); // mark the neigbours despite the direction
            residual[edge.where].push_back(node);
        }
    }

    int max_flow = 0, n_flow;

    while (_flowBfs(residual, source, dest, parent, capacity_matrix, flow_matrix) != -1) // are parinte
    {
      for (int k: residual[dest])
      {
          if (parent[k] != -1 && capacity_matrix[k][dest] > flow_matrix[k][dest])// facem flux doar cand stim ca muchia k dest nu e saturata
          {
               parent[dest] = k; // assume the parent is k
               n_flow = INT_MAX;

               int  n = dest;
               while (n != source)
               {
                    n_flow = min(n_flow, capacity_matrix[parent[n]][n] - flow_matrix[parent[n]][n] ); // make the flow of the path
                    n = parent[n];
               }

               if (n_flow > 0) // if it actually makes sense to update
               {
                   n = dest;
                   while (n != source)
                   {
                        flow_matrix[parent[n]][n] += n_flow; //update what the parent sends to child
                        flow_matrix[n][parent[n]] -= n_flow; // update what the child gets from parent
                        n = parent[n];
                   }

                   max_flow += n_flow;
               }
          }
      }
   }
    return max_flow;
}

vector<vector<int>> Graph::RoyFloyd(ifstream& in)
{
    vector<vector<int>> dist(_nr_vertex + 1, vector<int>(_nr_vertex + 1, INT_MAX));

    int aux;
    for (int node = 1; node <=_nr_vertex; ++node) // reading cost matrix
    {
        for (int i = 1; i <= _nr_vertex; ++i)
        {
           in >> aux;
           if (i != node && aux == 0)
                dist[node][i] = INT_MAX; // if there isnt a direct arc we assume the node is unreachable
           else dist[node][i] = aux;
        }
       dist[node][node] = 0;
    }
    for (int k = 1; k <= _nr_vertex; ++k)
    {
        for (int i = 1; i <= _nr_vertex; ++i)
        {
            for (int j = 1; j <= _nr_vertex; ++j)
            {
                if (dist[i][k] != INT_MAX && dist[k][j] != INT_MAX && dist[i][j] > dist[i][k] + dist[k][j])
                // update cost from i to j by calculating the costs from a  intermediate node
                    dist[i][j] = dist[i][k] + dist[k][j];
            }
        }
    }

    return dist;
}

int Graph::Darb(ifstream &in)
{
    int last, aux;
    _bfs(1, last);  // last = last node that gets max path length - so it's the furthest leaf from the root
    return  _bfs(last, aux); // from last we go back to the root and then descend into the other furthest leaf from root
}

vector<int> Graph::BFS(const int& S)
{
    /// Solving BFS from infoarena
    vector<int> result = _bfs(S);
    return vector<int> (result.begin() + 1, result.end());
}

int Graph::CompConexe()
{
     /// Solving DFS from infoarena
    int result = 0;
    vector<int> components(_nr_vertex + 1, -1);
    stack<int> aux; // not needed in the solving this

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (components[i] == -1)
        {
            result++;
            _dfs(i, result, components, aux);
        }
    }
    return result;
}

vector<int> Graph::sortareTopologica()
{
     /// Solving Sortare Topologica  from infoarena
    vector<int> components(_nr_vertex + 1, -1);
    vector<int> sol;
    stack<int> aux;

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (components[i] == -1) // A new component
        {
            _dfs(i, 0, components, aux); // going through DF tree
        }
    }

    while (!aux.empty())
    {
        sol.push_back(aux.top());
        aux.pop();
    }

    return sol; // vector of vertex's in sorted order
}

pair<int,vector<vector<int>>> Graph::compBiconexe()
{
    /// Solving Biconex  from infoarena
    vector<vector<int>> sol;
    vector<int> parent(_nr_vertex + 1, -1);
    vector<int> level(_nr_vertex + 1, -1);
    vector<int> rtr_lvl(_nr_vertex + 1, -1);
    stack<Edge> st;
    vector<vector<int>> crt_edg;

    _leveled_dfs(1, parent, level, rtr_lvl, st, sol, crt_edg);

    return pair<int, vector<vector<int>>> (sol.size(), sol); // return (numer of biconex components, sol[i] - vertex's of component)
}

vector<vector<int>> Graph::criticalConnections()
{
    /// Solving Critical Connections from leetcode
    // With a leveled DFS find critical edges on a conex non-oriented graph
    vector<vector<int>> sol;
    vector<int> parent;
    vector<int> level;
    vector<int> rtr_lvl;
    stack<Edge> st;
    vector<vector<int>> crt_edg;

    _leveled_dfs(0, parent, level, rtr_lvl, st, sol, crt_edg);
    return crt_edg; // crt_edg[i] - crt_edg[i][0] crt_edg[i][1] = critical edge
}

pair <int, vector<vector<int>>> Graph::componenteTareConexe()
{
    /// Solving CTC from infoarena
    vector<int> in_time(_nr_vertex + 1, -1);
    vector<int> time_return_vert(_nr_vertex +1, -1);
    stack<int> connection;
    vector<bool> in_connection(_nr_vertex + 1, 0);
    int time;
    time = 0;
    int sol_nr = 0;
    vector<vector<int>> sol;

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (in_time[i] == -1)
            _tarjan(i, time, in_time, time_return_vert, connection, in_connection, sol_nr, sol);
    }

    return make_pair(sol_nr, sol);
}

tuple<int, int, vector<Edge>> Graph::arborePartialCostMinim()
{
    pair<int, vector<Edge>> solution = Prim();
    return make_tuple(solution.first, solution.second.size(), solution.second);
}

void Graph::HavelHakimi(ifstream &in)
{
    /// Solving the problem "determine if array of degrees is a valid non-oriented graph" with Havel Hakimi theorem
    bool breakflag = 0;
    int aux, n, m = 0;
    int  nr_zeros = 0;
    vector<pair<int,int>> degrees;
    in >> n;
    for (int i = 1; i <= n; ++i)
    {
        in >> aux;
        if (aux > n) // if degree > number of vertex's => impossible
        {
            breakflag = 1;
            break;
        }
        if (!aux) nr_zeros++;
        m += aux;
        degrees.push_back(make_pair(i, aux));
    }
    if  ( m% 2 == 1 || breakflag ) cout << "Not a graph.\n";
    else   if (nr_zeros == n)
    {
        _resize(n, 0, 0, 0, 0);
        cout <<"An empty graph with " << n << " vertex's.\n";
    }
    else
     {
            bool answer = _HavelHakimi(n,m,degrees);
            if  (!answer) cout << "Not a graph.\n";
            else // valid graph so display edges
            for (int i = 1; i <= n; ++i)
            {
                cout << i <<": ";
                for (auto &e:_adjacency[i])
                {
                    cout << e.where << " ";
                }
                cout<< endl;
            }
     }
}

vector<bool> Graph::paduriDisjuncte(int N, int M, ifstream& in)
{
    vector<Forest> f = _makeForest(N);
    int x, y, cod;

    vector<bool> solution;
    while (M)
    {
        M--;

        in >> cod >> x >> y;
        if (cod == 1)   _unionForest(x, y, f);
        if (cod == 2)   solution.push_back(_checkForest(x,y,f));
    }

    return solution;
}

vector<int> Graph::BellmanFord()
{
    vector<int> rez = _BellmanFord(1);
    if (rez.size() != (unsigned int) _nr_vertex + 1)
        return  vector<int>(1,-1);
    else return  vector<int>(rez.begin()  + 2 , rez.end());
}

void Graph::_initializeAdMat(vector<vector<int>>& matrix)
{
    // Make copy of adjacency matrix
    matrix.resize(_nr_vertex + 1);

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        for (auto &edge: _adjacency[i])
            matrix[i].push_back(edge.where);
    }
}

vector<int> Graph::EulerianCycle()
{
    // Reunion of disjunct cycles
    //slaystudy.com/hierholzers-algorithm/
    vector<vector<int>> matrix;
    _initializeAdMat(matrix); //make another adjacency matrix not to destroy current one

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (matrix[i].size() % 2 == 1) // theorem
            return vector<int>(1, -1);
    }

    stack<int> nodes; // stack of nodes to make recursive alg. iterative
    nodes.push(1); // start from whatever node (bc cycle)

    vector<int> cycle; // solution
    int current_node, next_node, nr = 0;
    while (!nodes.empty()) // simulate recursion
    {
        current_node = nodes.top(); // take the current node

        if (nr == _nr_edge) // got eulerian cycle so break loop
            break;
        if (matrix[current_node].empty()) // passed through all neighbors
        {
            nr++; // mark nr of edges used
            nodes.pop();
            cycle.push_back(current_node); // add node to cycle
        }
        else
        {
                next_node = matrix[current_node].back(); // get a random neighbor
                matrix[current_node].pop_back(); // remove the edge from current node

                std::vector<int>::iterator position = find(matrix[next_node].begin(), matrix[next_node].end(), current_node);
                if (position != matrix[next_node].end())
                    matrix[next_node].erase(position); // remove the edge from the other node side
                nodes.push(next_node); // add node to recursion
        }
    }

    if (nr != _nr_edge) // there's no cycle.
        return vector<int> (1, -1);

    return cycle;
}

void infoarenaBFS()
{
    ifstream in("bfs.in");
    int n, m, s;
    in >> n >> m >> s;
    Graph g(n,m, 1);
    g.readEdges(in);
    vector<int> sol = g.BFS(s);

    ofstream out("bfs.out");
    for (unsigned int i = 0; i < sol.size(); ++i)
        out << sol[i] << " ";
    out.close();
}

void infoarenaDFS()
{
    ifstream in("dfs.in");
    int n, m;
    in >> n >> m;
    Graph g(n,m, 0);
    g.readEdges(in);
    in.close();
    int sol = g.CompConexe();

    ofstream out("dfs.out");
    out << sol;
    out.close();
}

void infoarenaSortareTopologica()
{
    ifstream in("sortaret.in");
    ofstream out("sortaret.out");
    int n, m;
    in >> n >> m;
    Graph g(n,m,1);
    g.readEdges(in);
    in.close();

    vector<int> sol = g.sortareTopologica();
    for (unsigned int i = 0; i < sol.size(); ++i)
        out << sol[i] << " ";
    out.close();
}

void infoarenaBiconex()
{
     ifstream in("biconex.in");
    ofstream out("biconex.out");
    int n, m;
    in >> n >> m;
    Graph g(n,m,0);
    g.readEdges(in);
    in.close();
    pair<int, vector<vector<int>>> sol = g.compBiconexe();
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
}

void leetCriticalConnections()
{
     ifstream in("criticalconnections.in");
    ofstream out("criticalconnections.out");
    int n, m;
    in >> n >> m;
    Graph g(n,m,0);
    g.readEdges(in);
    in.close();
    vector<vector<int>> crt_edg = g.criticalConnections();

    for (unsigned int i = 0; i < crt_edg.size(); ++i)
    {
        for (unsigned int e_idx = 0; e_idx < crt_edg[i].size(); ++e_idx)
        {
            out << crt_edg[i][e_idx] << " ";
        }
        out << "\n";
    }
}

void infoarenaCTC()
{
    ifstream in("ctc.in");
    int n, m;
    in >> n >> m;
    Graph g(n,m,1);
    g.readEdges(in);
    in.close();
    pair<int, vector<vector<int>>> solution = g.componenteTareConexe();
    ofstream out("ctc.out");
    out << solution.first << endl;
    for (int i = 0; i < solution.first; ++i)
    {
        for (auto &it: solution.second[i])
            out << it << " ";
        out << "\n";
    }
    out.close();

}

void solveHH()
{
    ifstream in("havelhakimi.in");
    Graph g(0,0, 0);
    g.HavelHakimi(in);
}

void infoarenaAPM()
{
    ifstream in("apm.in");
     int n, m;
    in >> n >> m;
    Graph g(n,m,0,1);
    g.readEdges(in);
    in.close();
    ofstream out("apm.out");
    tuple<int,int, vector<Edge>> sol = g.arborePartialCostMinim();
    out << get<0>(sol) << "\n" << get<1>(sol) << "\n";

    for (auto &edge:get<2>(sol))
    {
        out << edge.from << " " << edge.where << "\n";
    }

    out.close();
}

void  infoarenaDijkstra()
{
    ifstream in("dijkstra.in");
    int n, m;
    in >> n >> m;
    Graph g(n,m,1,1);
    g.readEdges(in);
    in.close();
    vector<int> solution = g.Dijkstra(1);
    ofstream out("dijkstra.out");
    for (int i = 0;  i < n - 1 ; ++i)
    {
        if (solution[i] == INT_MAX) out << 0 << " ";
        else out << solution[i] << " ";
    }
    out.close();
}

void infoarenaDisjunct()
{
    ifstream in("disjoint.in");
    int n, m;
    in >> n >> m;
    Graph g;
    vector<bool> solution = g.paduriDisjuncte(n, m, in);
    in.close();
    ofstream out("disjoint.out");
    for (auto it:solution) it ? out <<"DA\n" : out << "NU\n";
    out.close();
}

void infoarenaBellmanFord()
{
     ifstream in("bellmanford.in");
    int n, m;
    in >> n>> m;
    Graph g(n,m,1,1);
    g.readEdges(in);
    in.close();
    ofstream out("bellmanford.out");
    vector<int> sol = g.BellmanFord();
    if (sol.size() != (unsigned int) n - 1)
        out << "Ciclu negativ!\n";
    else
    {
        for (auto &it:sol)
            out << it << " ";
    }
    out.close();
}

void infoarenaMaxflow()
{
    ifstream in("maxflow.in");
    int n, m;
    in >> n >> m;
    Graph g(n, m, 1, 0, 1);
    g.readEdges(in);
    in.close();
    ofstream out("maxflow.out");
    out << g.maxFlow(1, n);
    out.close();
}

void infoarenaRoyFloyd()
{
     ifstream in("royfloyd.in");
    int n;
    in >> n;
    Graph g(n, n*n, 1, 1, 0);
    vector<vector<int>> sol =  g.RoyFloyd(in);
    in.close();
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
}

void infoarenaDarp()
{
     ifstream in("darb.in");
    ofstream out("darb.out");
    int n;
    in >> n;
    Graph g(n, n - 1, 0, 0, 0);
    g.readEdges(in);
    out << g.Darb(in);
    in.close();
    out.close();
}

int main()
{
    ifstream in("ciclueuler.in");
    ofstream out("ciclueuler.out");

    int n, m;
    in >> n >> m;
    Graph g(n, m);
    g.readEdges(in);
    vector<int> sol = g.EulerianCycle();
    if (sol.size() == 1 && sol[0] == -1)
        out << "-1\n";
    else
    {
        for (auto &node: sol)
            out << node << " ";
    }
    return 0;
}
