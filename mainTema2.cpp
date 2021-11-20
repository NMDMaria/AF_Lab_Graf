#include <bits/stdc++.h>

using namespace std;

struct Edge
{
    int from, where;
    int cost;
};

struct Weighted_edge
{
    int where, weight;
    bool operator()(const Weighted_edge &left, const Weighted_edge& right)
    {
        return left.weight > right.weight;
    }
};

bool pairsortsnddesc(const pair<int,int>& i, const pair<int,int>& j)
{
    return i.second > j.second;
}

class Graph
{
 private:
    int _nr_vertex;
    int _nr_edge;
    bool _oriented;
    bool _weighted;

    vector<vector<Weighted_edge>> _adjacency;
    // Internal methods.
    ifstream& _readEdges(ifstream&);
    vector<int> _popEdges(stack<Edge>& , const Edge&);
    void _popVertex(const int&, stack<int>&, int&, vector<vector<int>>&, vector<bool>&);
    vector<int> _bfs(const int&);
    void _dfs(const int&, const int&, vector<int>&, stack<int>&);
    void _leveled_dfs(const int&, vector<int>& , vector<int>&, vector<int>& , stack<Edge>&, vector<vector<int>>&, vector<vector<int>>&);
    void _tarjan(const int&, int&, vector<int>&, vector<int>&, stack<int>&, vector<bool>& , int&, vector<vector<int>>&);
    void _HavelHakimi(const int&, const int&, vector<pair<int,int>>& , bool&);
    void _addEdge(const Edge&);
    void _resize(const int&, const int&, const bool&, const bool&);
    pair<int,vector<Edge>> _Prim();

 public:
    // Constructors
    Graph(const int&, const int&, const bool&, const bool&, const vector<Edge>&);
    Graph(const int& n = 0, const int& m = 0, const bool& o = 0, const bool& w = 0);

    // Functions and methods to solve the requirements
    vector<int> solveBFS(ifstream&);
    int solveDFS(ifstream&);
    vector<int> solveTopo(ifstream&);
    pair<int,vector<vector<int>>> solveBiconex(ifstream&);
    vector<vector<int>> criticalConnections(ifstream&);
    pair<int, vector<vector<int>>> solveCTC(ifstream &);
    void HavelHakimi(ifstream &);
    void solveAPM(ifstream &);
};

/// Constructors
Graph::Graph(const int& n, const int& m, const bool& o, const bool& w, const vector<Edge>& edges)
{


   _nr_vertex= n;
   _nr_edge = m;
   _oriented = o;
   _weighted = w;
   if (n != 0)
   {
       _adjacency.resize(n + 1);
       for (auto &e: edges)
       {
           _adjacency[e.from].push_back({e.where, e.cost});
           if (!_oriented) _adjacency[e.where].push_back({e.from, e.cost});
       }
   }
}

Graph::Graph(const int& n, const int& m, const bool& o, const bool& w)
{

   _nr_vertex= n;
   _nr_edge = m;
   _oriented = o;
   _weighted = w;
   if (n != 0) _adjacency.resize(n + 1);
}

/// Definitions of internal methods
void Graph::_resize(const int& new_n, const int& new_m, const bool& new_o, const bool& new_w)
{

    _nr_vertex= new_n;
    _nr_edge = new_m;
    _oriented = new_o;
    _weighted = new_w;
    _adjacency.clear();
    _adjacency.resize(new_n + 1);
}

void Graph::_addEdge(const Edge& e)
{
    /// Add an edge to the graph
    _adjacency[e.from].push_back({e.where, e.cost});
    if (!_oriented) _adjacency[e.where].push_back({e.from, e.cost});
}

ifstream& Graph::_readEdges(ifstream& in)
{
    /// Reads the edges of the graph
    // in = pointer of the location of the edges in the open for reading file

    int x, y, c;
    for(int i = 0; i < _nr_edge; ++i)
    {
        in >> x >> y;
        if (_weighted)
        {
            in >> c;
            _addEdge({x, y, c}); // Creates pair of type Edge
        }
        else _addEdge({x, y, 0});
    }
    return in;
}

vector<int> Graph::_bfs(const int& startV)
{
    /// Solves BFS  returning the distances vector
    // startV = the start vertex of the bfs
    // return: distances = vector of _nr_vertex + 1 length with the distances from the start vertex to the i-th one

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

void Graph::_dfs(const int& start, const int& marker, vector<int>& mark, stack<int>& exit_time)
{
    /// Solves DFS recursively
    // exit_time = stack in the order of exiting vertex's DFS

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
    // sol_nr = number of strong connected components
    // sol[i] = vertex's of the strong connected component
    // in[i] = the vertex i was visited in the current component

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
    // start = current vertex in the dfs
    // level [i] = depth of the vertex in the dfs tree
    // parent[i] = the parent of the vertex i in the dfs tree
    // return_level[i] = the level that the vertex i can return to using return edges
    // expl_edges = stack of edges in the order of discovery - used to discover the biconex components
    // biconex_comps = vector of  the biconex components and their vertex's
     // critical_edges[i] = vector of 2 elements (leetcode restriction) that signify a critical edge

    if ( parent.size() != (unsigned) (_nr_vertex + 1) && level.size() !=(unsigned) (_nr_vertex + 1)  &&  return_level.size() !=(unsigned) (_nr_vertex + 1) )
    {
        //  First iteration of function, initialization of vectors
        parent.resize(_nr_vertex + 1);
        parent.assign(_nr_vertex + 1, -1);
        level.resize(_nr_vertex + 1);
        return_level.resize(_nr_vertex + 1);
        level.assign(_nr_vertex + 1, -1);
        return_level.assign(_nr_vertex + 1, -1);
        parent[start] = 0; // The root of the DFS tree
        level[start] = 0;
        return_level[start] = 0;
    }

    int nr_children = 0;
    for (auto &child: _adjacency[start])
    {
        nr_children ++; // mark as child
        if (parent[child.where] == -1)
        {
            expl_edges.push({start, child.where, 0});
            parent[child.where] = start;
            level[child.where] = level[start] + 1;
            return_level[child.where] = level[child.where]; // now the lowest level  reached  from the child is his exact level

            _leveled_dfs(child.where, parent, level, return_level, expl_edges, biconex_comps, critical_edges);

            return_level[start] = min(return_level[start], return_level[child.where]); // passed the rest of the DF tree of the child
            // can modify the lowest level reached in case there was a vertex with a return edge

            if (return_level[child.where] > level[start]) // the child cant return higher then his parent
            {
                critical_edges.push_back(vector<int>{start,child.where, 0});
            }

            if (parent[start] == 0 && nr_children >= 2)
            {
                vector<int> new_comp = _popEdges(expl_edges, {start,child.where});  // get the biconex component of start
                biconex_comps.push_back(new_comp);
            }

            if (parent[start] != 0 && return_level[child.where] >= level[start]) // the child can return to a lower level then his parent
            {
                vector<int> new_comp = _popEdges(expl_edges, {start,child.where}); // get th ebiconex component of start
                biconex_comps.push_back(new_comp);
            }
        }

        else if (child.where != parent[start]) return_level[start] = min(return_level[start], level[child.where]); // update the lowest level reachable with return edges
    }
}

void Graph::_tarjan(const int& start, int& time, vector<int>& in_time, vector<int>& time_return_vert, stack<int>& connection, vector<bool>& in_connection, int& nr_ctc, vector<vector<int>>& strong_connected)
{
    /// Given a vertex determine his strong connected component recursively (DF based)
    // start = current vertex
    // time = the nr. of  iterations, meaning the time of discovery of a vertex in DF tree
    // in_time[i] = discovery time of i
    // time_return_vert[i] = the lowest discovery time in the i DF sub-tree  + return edges
    // in_connection[i] = true if i is in the current SCC
    // strong_connected & nr_ctc = nr of strong connected components and the respective vertex's

    time++; // mark a new iteration
    in_time[start] = time;
    time_return_vert[start] = time; // no return edge known
    in_connection[start] = 1;
    connection.push(start);

    for (auto &child: _adjacency[start])
    {
        if (in_time[child.where] == -1)
        {
            _tarjan(child.where, time,  in_time, time_return_vert, connection, in_connection, nr_ctc, strong_connected);
            time_return_vert[start] = min(time_return_vert[start], time_return_vert[child.where]); // update in case a return edge was found in the child's  df sub-tree
        }
        else if(in_connection[child.where]) time_return_vert[start] = min(time_return_vert[start], in_time[child.where]);
    }

    if (time_return_vert[start] == in_time[start]) // a vertex that doesn't have a return edge => end of SCC
    {
        _popVertex(start, connection, nr_ctc, strong_connected, in_connection);
    }
}

void Graph::_HavelHakimi(const int& n, const int& nr_d, vector<pair<int,int>>& degrees , bool& breakflag)
{
    /// Given a vector of pairs (Vertex, degree) determine if it is a valid graph
    int  m = nr_d/2; // number of edges is sum of degrees/2
    _resize(n, m, 0, 0);
    breakflag = 0;
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
                if (degrees[k].second < 0) { breakflag = 1; break;}
                m--;

                _addEdge({max_dg.first, degrees[k].first});
            }
        }
        if (max_dg.second != 0) {breakflag = 1; break;}
    }

    if (m!=0) breakflag = 1; // Couldn't add the number of desired edges so not a valid graph
}

pair<int, vector<Edge>> Graph::_Prim()
{
    priority_queue<Weighted_edge, vector<Weighted_edge>, Weighted_edge> heap;
    vector<int> cost(_nr_vertex + 1, INT_MAX);
    vector<int> parent(_nr_vertex + 1, -1);
    vector<bool> in_APM(_nr_vertex + 1, 0);
    vector<Edge> edges;

    int node, total_cost = 0;
    cost[1] = 0;

    heap.push({1, cost[1]});

    while (!heap.empty())
    {
        node = heap.top().where;
        heap.pop();

        if (!in_APM[node])
        {
            in_APM[node] = 1;
            total_cost += cost[node];

            for (auto &neighbour: _adjacency[node])
            {
                if (cost[neighbour.where] > neighbour.weight && !in_APM[neighbour.where])
                {
                    heap.push({neighbour.where, neighbour.weight});
                    parent[neighbour.where] = node;
                    cost[neighbour.where] = neighbour.weight;
                }
            }

            if (parent[node] != -1)  edges.push_back({parent[node], node, cost[node]});
        }
    }

    return make_pair(total_cost, edges);
}

/// Procedures for solving the requirements
vector<int> Graph::solveBFS(ifstream &in)
{
    /// Solving BFS from infoarena
    int s;
    in >> s;
    _readEdges(in);
    in.close();
    vector<int> result = _bfs(s);
    return vector<int> (result.begin() + 1, result.end());
}

int Graph::solveDFS(ifstream& in)
{
     /// Solving DFS from infoarena
    _readEdges(in);
    in.close();
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

vector<int> Graph::solveTopo(ifstream& in)
{
     /// Solving Sortare Topologica  from infoarena
    _readEdges(in);
    in.close();
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

pair<int,vector<vector<int>>> Graph::solveBiconex(ifstream &in)
{
    /// Solving Biconex  from infoarena
     _readEdges(in);
    in.close();
    vector<vector<int>> sol;
    vector<int> parent;
    vector<int> level;
    vector<int> rtr_lvl;
    stack<Edge> st;
    vector<vector<int>> crt_edg;

    _leveled_dfs(1, parent, level, rtr_lvl, st, sol, crt_edg);

    vector<int> last_cmp;
    int x,y;
    while (!st.empty()) // Getting the last biconex component that remains in the stack
    {
        x = st.top().from;
        y = st.top().where;
        st.pop();
        last_cmp.push_back(x);
        last_cmp.push_back(y);
    }
    sol.push_back(last_cmp);

    return pair<int, vector<vector<int>>> (sol.size(), sol); // return (numer of biconex components, sol[i] - vertex's of component)
}

vector<vector<int>> Graph::criticalConnections(ifstream &in)
{
    /// Solving Critical Connections from leetcode
    // With a leveled DFS find critical edges on a conex non-oriented graph
    _readEdges(in);
    in.close();
    vector<vector<int>> sol;
    vector<int> parent;
    vector<int> level;
    vector<int> rtr_lvl;
    stack<Edge> st;
    vector<vector<int>> crt_edg;

    _leveled_dfs(0, parent, level, rtr_lvl, st, sol, crt_edg);
    return crt_edg; // crt_edg[i] - crt_edg[i][0] crt_edg[i][1] = critical edge
}

pair <int, vector<vector<int>>> Graph::solveCTC(ifstream &in)
{
    /// Solving CTC from infoarena
     _readEdges(in);
    in.close();
    vector<int> in_time;
    vector<int> time_return_vert;
    stack<int> connection;
    vector<bool> in_connection;
    in_connection.resize(_nr_vertex + 1);
    in_connection.assign(_nr_vertex + 1, 0);
    int time;
    in_time.resize(_nr_vertex + 1);
    in_time.assign(_nr_vertex + 1, -1);
    time_return_vert.resize(_nr_vertex + 1);
    time_return_vert.assign(_nr_vertex + 1, -1);
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

void Graph::solveAPM(ifstream &in)
{
    int n, m;
    in >> n >> m;
    _resize(n, m, 0, 1);
    _readEdges(in);
    pair<int, vector<Edge>> solution = _Prim();
    ofstream out("apm.out");
    out << solution.first << "\n" << n - 1 << "\n";

    for (unsigned int i = 0;  i < solution.second.size(); ++i)
        out << solution.second[i].from << " " << solution.second[i].where << "\n";

    out.close();
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
        _resize(n, 0, 0, 0);
        cout <<"An empty graph with " << n << " vertex's.\n";
    }
    else
     {
            _HavelHakimi(n,m,degrees, breakflag);

            if  (breakflag) cout << "Not a graph.\n";
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


void infoarenaBFS()
{
    ifstream in("bfs.in");
    int n, m;
    in >> n >> m;
    Graph g(n,m, 1);
    vector<int> sol = g.solveBFS(in);

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
    int sol = g.solveDFS(in);

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
    vector<int> sol = g.solveTopo(in);
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
    pair<int, vector<vector<int>>> sol = g.solveBiconex(in);
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
    vector<vector<int>> crt_edg = g.criticalConnections(in);

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
    pair<int, vector<vector<int>>> solution = g.solveCTC(in);
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

void solveAPM()
{
    ifstream in("apm.in");
    Graph g;
    g.solveAPM(in);
}

int main()
{
    solveAPM();
    return 0;
}
