#include <bits/stdc++.h>

using namespace std;

struct Edge
{
    /// Struct that helps with the definition of pairs for the edges of the graphs
    int from, where;
};

bool pairsortsnddesc(const pair<int,int>& i, const pair<int,int>& j)
{
    /// Function used to sort a vector of pairs descending by the second component
    return i.second > j.second;
}

class Graph
{
public:
    Graph(const int&, const int&, const bool&, const vector<Edge>&);
    Graph(const int& n = 0, const int& m = 0, const bool& o = 0);

    vector<int> SolveBFS(ifstream&);
    int SolveDFS(ifstream&);
    vector<int> SolveTopo(ifstream&);
    pair<int,vector<vector<int>>> SolveBiconex(ifstream&);
    vector<vector<int>> CriticalConnections(ifstream&);
    pair<int, vector<vector<int>>> SolveCTC(ifstream &);
    void HavelHakimi(ifstream &);

 private:
    int _nr_vertex;
    int _nr_edge;
    bool _oriented;
    vector<vector<int>> _adjacency;

    ifstream& _readEdges(ifstream&);
    vector<int> _popEdges(stack<Edge>& , const Edge&);
    vector<int>  _popVertex(const int&, stack<int>&,  vector<bool>&);
    vector<int> _bfs(const int&);
    void _dfs(const int&, const int&, vector<int>&, stack<int>&);
    void _leveled_dfs(const int&, vector<int>& , vector<int>&, vector<int>& , stack<Edge>&, vector<vector<int>>&, vector<vector<int>>&);
    void _tarjan(const int&, int&, vector<int>&, vector<int>&, stack<int>&, vector<bool>& , int&, vector<vector<int>>&);
    bool _havel_hakimi(const int&, const int&, vector<pair<int,int>>& );
    void _addEdge(const Edge&);
    void _resize(const int&, const int&, const bool&);
};

Graph::Graph(const int& n, const int& m, const bool& o, const vector<Edge>& edges)
{
   _nr_vertex= n;
   _nr_edge = m;
   _oriented = o;
   if (n != 0)
   {
       _adjacency.resize(n + 1);
       for (auto &e: edges)
       {
           _adjacency[e.from].push_back(e.where);
           if (!_oriented) _adjacency[e.where].push_back(e.from);
       }
   }
}

Graph::Graph(const int& n, const int& m, const bool& o)
{
   _nr_vertex= n;
   _nr_edge = m;
   _oriented = o;
   if (n != 0)   _adjacency.resize(n + 1);
}

void Graph::_resize(const int& new_n, const int& new_m, const bool& new_o)
{
    _nr_vertex= new_n;
    _nr_edge = new_m;
    _oriented = new_o;
    _adjacency.clear();
    _adjacency.resize(new_n + 1);
}

void Graph::_addEdge(const Edge& e)
{
    _adjacency[e.from].push_back(e.where);
    if (!_oriented) _adjacency[e.where].push_back(e.from);
}

ifstream& Graph::_readEdges(ifstream& in)
{
    int x, y;
    for(int i = 0; i < _nr_edge; ++i)
    {
        in >> x >> y;
        _addEdge({x,y}); // Creates pair of type Edge
    }
    return in;
}

vector<int> Graph::_bfs(const int& startv)
{
    /// Solves BFS  returning the distances vector
    vector<int> distances( _nr_vertex + 1, -1);
    distances[startv] = 0;
    queue<int> order_list;
    order_list.push(startv);
    int current;

    while (!order_list.empty())
    {
        current = order_list.front();
        for (auto &neighbour: _adjacency[current])
        {
            if (distances[neighbour] == -1)
            {
                distances[neighbour] = distances[current] + 1;
                order_list.push(neighbour);
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
        if (mark[neighbour] == -1)
            _dfs(neighbour, marker, mark, exit_time);
    }
    exit_time.push(start);
}

vector<int> Graph::_popEdges(stack<Edge>& st, const Edge& last_edg)
{
    /// Given an stack of edges pop the elements of the stack till we reach the one given
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

vector<int> Graph::_popVertex(const int& start, stack<int>& st, vector<bool> &in)
{
    /// Given a stack of vertex's  pop the elements of the stack till we reach the one given
    int aux;
    vector<int> aux2;
    do
    {
        aux = st.top();
        in[aux] = 0;
        st.pop();
        aux2.push_back(aux);
    }while (aux != start);
    return aux2;
}

void Graph::_leveled_dfs(const int& start, vector<int>& parent, vector<int>& level, vector<int>&  return_level, stack<Edge>& expl_edges, vector<vector<int>>& biconex_comps, vector<vector<int>>& critical_edges)
{
    /// Given a start vertex do a recursive modified DFS, with the purpose of finding the critical vertex's/edges (similar to Tarjan's algorithm)
    if ( parent.size() != (unsigned) (_nr_vertex + 1) && level.size() !=(unsigned) (_nr_vertex + 1)  &&  return_level.size() !=(unsigned) (_nr_vertex + 1) )
    {
        parent.resize(_nr_vertex + 1);
        parent.assign(_nr_vertex + 1, -1);
        level.resize(_nr_vertex + 1);
        return_level.resize(_nr_vertex + 1);
        level.assign(_nr_vertex + 1, -1);
        return_level.assign(_nr_vertex + 1, -1);
        parent[start] = 0;
        level[start] = 0;
        return_level[start] = 0;
    }

    int nr_children = 0;
    for (auto &child: _adjacency[start])
    {
        nr_children ++;
        if (parent[child] == -1)
        {
            expl_edges.push({start, child});
            parent[child] = start;
            level[child] = level[start] + 1;
            return_level[child] = level[child];

            _leveled_dfs(child, parent, level, return_level, expl_edges, biconex_comps, critical_edges);

            return_level[start] = min(return_level[start], return_level[child]); // passed the rest of the DF tree of the child
            // can modify the lowest level reached in case there was a vertex with a return edge

            if (return_level[child] > level[start]) // the child cant return higher then his parent
            {
                critical_edges.push_back(vector<int>{start,child}); // so this edge is critical in the graph
            }

            if (parent[start] == 0 && nr_children >= 2)
            {
                vector<int> new_comp = _popEdges(expl_edges, {start,child});  // get the biconex component of start
                biconex_comps.push_back(new_comp);
            }

            if (parent[start] != 0 && return_level[child] >= level[start])
            {
                vector<int> new_comp = _popEdges(expl_edges, {start,child}); // get the biconex component of start
                biconex_comps.push_back(new_comp);
            }
        }

        else if (child != parent[start]) return_level[start] = min(return_level[start], level[child]); // update the lowest level reachable with return edges
    }
}

void Graph::_tarjan(const int& start, int& time, vector<int>& in_time, vector<int>& time_return_vert, stack<int>& connection, vector<bool>& in_connection, int& nr_ctc, vector<vector<int>>& strong_connected)
{
    /// Given a vertex determine his strong connected component recursively (DF)
    time++;
    in_time[start] = time;
    time_return_vert[start] = time;
    in_connection[start] = 1;
    connection.push(start);

    for (auto &child: _adjacency[start])
    {
        if (in_time[child] == -1)
        {
            // continue df
            _tarjan(child, time,  in_time, time_return_vert, connection, in_connection, nr_ctc, strong_connected);
            time_return_vert[start] = min(time_return_vert[start], time_return_vert[child]); // update in case a return edge was found in the child's  df sub-tree
        }
        else if(in_connection[child]) time_return_vert[start] = min(time_return_vert[start], in_time[child]);
    }

    if (time_return_vert[start] == in_time[start]) // a vertex that doesn't have a return edge => end of SCC
    {
        vector<int> aux;
        aux = _popVertex(start, connection, in_connection);
        nr_ctc++;
        strong_connected.push_back(aux);
    }
}

bool Graph::_havel_hakimi(const int& n, const int& nr_d, vector<pair<int,int>>& degrees )
{
    /// Given a vector of pairs (Vertex, degree) determine if it is a valid graph
    int  m = nr_d/2;
    _resize(n, m, 0);
    vector<vector<bool>> matrix_adj;
    matrix_adj.resize(n + 1);
    for (int i = 1; i <= n; ++i)
        matrix_adj[i].resize(n + 1, 0);
    while (m)
    {
        sort(degrees.begin(), degrees.end(), pairsortsnddesc);
        if (degrees.size() == 0 || degrees[0].second == 0) break;
        pair<int,int> max_dg;
        max_dg.first = degrees[0].first;
        max_dg.second = degrees[0].second;
        degrees.erase(degrees.begin());
        for (unsigned int k = 0; k < degrees.size() && max_dg.second > 0 ; ++k)
        {
            if (!matrix_adj[max_dg.first][degrees[k].first])
            {
                matrix_adj[max_dg.first][degrees[k].first] = matrix_adj[degrees[k].first][max_dg.first] = 1; // mark new edge
                degrees[k].second--; // lower degrees
                max_dg.second--;
                if (degrees[k].second < 0) return 0;
                m--;

                _addEdge({max_dg.first, degrees[k].first});
            }
        }
        if (max_dg.second != 0) return 0;
    }

    if (m!=0) return 0;
    return 1;
}

vector<int> Graph::SolveBFS(ifstream &in)
{
    /// Solving BFS from infoarena
    int s;
    in >> s;
    _readEdges(in);
    in.close();
    vector<int> result = _bfs(s);
    return vector<int> (result.begin() + 1, result.end());
}

int Graph::SolveDFS(ifstream& in)
{
     /// Solving DFS from infoarena
    _readEdges(in);
    in.close();
    int result = 0;
    vector<int> components(_nr_vertex + 1, -1);
    stack<int> aux;

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (components[i] == -1) // A new component
        {
            result++;
            _dfs(i, result, components, aux);
        }
    }
    return result;
}

vector<int> Graph::SolveTopo(ifstream& in)
{
     /// Solving Sortare Topologica  from infoarena
    _readEdges(in);
    in.close();
    vector<int> components(_nr_vertex + 1, -1);
    vector<int> sol;
    stack<int> aux;

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (components[i] == -1)
        {
            _dfs(i, 0, components, aux);
        }
    }

    while (!aux.empty())
    {
        sol.push_back(aux.top());
        aux.pop();
    }

    return sol; // vector of vertex's in sorted order
}

pair<int,vector<vector<int>>> Graph::SolveBiconex(ifstream &in)
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

    _leveled_dfs(1, parent, level, rtr_lvl, st, sol, crt_edg); // Graph is connex so we can start from anywhere

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

vector<vector<int>> Graph::CriticalConnections(ifstream &in)
{
    /// Solving Critical Connections from leetcode
    _readEdges(in);
    in.close();
    vector<vector<int>> sol;
    vector<int> parent;
    vector<int> level;
    vector<int> rtr_lvl;
    stack<Edge> st;
    vector<vector<int>> crt_edg;

    _leveled_dfs(0, parent, level, rtr_lvl, st, sol, crt_edg); // Graph is conex start from anywhere
    return crt_edg; // crt_edg[i] - crt_edg[i][0] crt_edg[i][1] = critical edge
}

pair <int, vector<vector<int>>> Graph::SolveCTC(ifstream &in)
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
    in_time.resize(_nr_vertex + 1);
    in_time.assign(_nr_vertex + 1, -1);
    time_return_vert.resize(_nr_vertex + 1);
    time_return_vert.assign(_nr_vertex + 1, -1);
    int time = 0, sol_nr = 0;
    vector<vector<int>> sol;

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (in_time[i] == -1)
            _tarjan(i, time, in_time, time_return_vert, connection, in_connection, sol_nr, sol);
    }

    return make_pair(sol_nr, sol);
}

void Graph::HavelHakimi(ifstream &in)
{
    /// Solving the problem "determine if array of degrees is a valid non-oriented graph" with Havel Hakimi theorem
    bool breakflag = 0;
    int aux, n, m = 0;
    int  nr_zeros = 0;
    vector<pair<int,int>> degrees;
    in >> n; // read number of vertex's
    for (int i = 1; i <= n; ++i)
    {
        in >> aux;
        if (aux > n) // if degree > number of vertex's => impossible
        {
            breakflag = 1;
            break;
        }
        if (!aux) nr_zeros++;
        m += aux; // the number of edges is the number of degrees/2
        degrees.push_back(make_pair(i, aux));
    }
    if  ( m % 2 == 1 || breakflag ) cout << "Not a graph.\n";
    else   if (nr_zeros == n) { _resize(n, 0, 0);  cout <<"An empty graph with " << n << " vertex's.\n"; }
    else
     {
            breakflag = _havel_hakimi(n,m,degrees);
            if  (!breakflag) cout << "Not a graph.\n";
            else // valid graph so display edges
            for (int i = 1; i <= n; ++i)
            {
                cout << i <<": ";
                for (auto &e:_adjacency[i])
                {
                    cout << e << " ";
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
    vector<int> sol = g.SolveBFS(in);

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
    int sol = g.SolveDFS(in);

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
    vector<int> sol = g.SolveTopo(in);
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
    pair<int, vector<vector<int>>> sol = g.SolveBiconex(in);
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
    vector<vector<int>> crt_edg = g.CriticalConnections(in);

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
    pair<int, vector<vector<int>>> solution = g.SolveCTC(in);
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

int main()
{
    // call a function from above
    solveHH();
    return 0;
}
