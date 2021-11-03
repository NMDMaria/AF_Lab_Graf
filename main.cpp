#include <bits/stdc++.h>

using namespace std;

struct Edge
{
    /// Struct that helps with the definition of pairs for the edges of the graphs
    // from = the start vertex of the edge
    // where = the end vertex of the edge
    int from, where;
};

class Graph
{
 private:
    int _nr_vertex;
    int _nr_edge;
    bool _oriented;
    vector<vector<int>> _adjacency;

    // Internal methods.
    ifstream& _readEdges(ifstream&);
    vector<int> _popEdges(stack<Edge>& , const Edge&);
    vector<int> _bfs(const int&);
    void _dfs(const int&, const int&, vector<int>&, stack<int>&);
    void _leveled_dfs(const int&, vector<int>& , vector<int>&, vector<int>& , stack<Edge>&, vector<vector<int>>&);
    void _addEdge(const Edge&);
    void _resize(const int&, const int&, const bool&, const vector<Edge>& edges);
    void _resize(const int&, const int&, const bool&);

 public:
    // Constructors
    Graph(const int&, const int&, const bool&, const vector<Edge>&);
    Graph(const int& n = 0, const int& m = 0, const bool& o = 0);

    // Functions and methods to solve the requirements
    vector<int> solveBFS(ifstream&);
    int solveDFS(ifstream&);
    vector<int> solveTopo(ifstream&);
    pair<int,vector<vector<int>>> solveBiconex(ifstream&);
};

/// Constructors
Graph::Graph(const int& n, const int& m, const bool& o, const vector<Edge>& edges)
{
    /// Constructor with full parameters = creates a graph with his edges
    // n = number of vertexes
    // m = number of edges
    // o = bool type that when true means the graph created is oriented
    // edges = vector of the edges in the graph retain by <from what vertex>, <to what vertex>

   _nr_vertex= n;
   _nr_edge = m;
   _oriented = o;
   if (n != 0) // no need to resize adjacency vector for a graph with no vertex
   {
       _adjacency.resize(n + 1);
       for (auto &e: edges)
       {
           _adjacency[e.from].push_back(e.where);
           if (!_oriented) _adjacency[e.where].push_back(e.from); // if the graph is not oriented can put the "reverse" edge to get access on both sides
       }
   }
}

Graph::Graph(const int& n, const int& m, const bool& o)
{
     /// Constructor with the basic parameters = creates a graph without the adjacency vector
    // n = number of vertexes
    // m = number of edges
    // o = bool type that when true means the graph created is oriented

   _nr_vertex= n;
   _nr_edge = m;
   _oriented = o;
   if (n != 0)   _adjacency.resize(n + 1);
}

/// Definitions of internal methods
void Graph::_resize(const int& new_n, const int& new_m, const bool& new_o, const vector<Edge>& new_edges)
{
    /// Transform the graph to match the description
    // See full parameter constructor for parameters description.
    _nr_vertex= new_n;
    _nr_edge = new_m;
    _oriented = new_o;
    _adjacency.clear();
     _adjacency.resize(new_n + 1);
    for (auto &e: new_edges)
   {
       _adjacency[e.from].push_back(e.where);
       if (!_oriented) _adjacency[e.where].push_back(e.from);
   }
}

void Graph::_resize(const int& new_n, const int& new_m, const bool& new_o)
{
    /// Transform the graph to match the description
    // See full parameter constructor for parameters description.
    _nr_vertex= new_n;
    _nr_edge = new_m;
    _oriented = new_o;
    _adjacency.clear();
    _adjacency.resize(new_n + 1);
}

void Graph::_addEdge(const Edge& e)
{
    /// Add an edge to the graph
    _adjacency[e.from].push_back(e.where);
    if (!_oriented) _adjacency[e.where].push_back(e.from);
}

ifstream& Graph::_readEdges(ifstream& in)
{
    /// Reads the edges of the graph
    // in = pointer of the location of the edges in the open for reading file

    int x, y;
    for(int i = 0; i < _nr_edge; ++i)
    {
        in >> x >> y;
        _addEdge({x,y}); // Creates pair of type Edge
    }
    return in;
}

vector<int> Graph::_bfs(const int& startV)
{
    /// Solves BFS  returning the distances vector
    // startV = the start vertex of the bfs
    // return: distances = vector of _nr_vertex + 1 length with the distances from the start vertex to the i-th one

    vector<int> distances( _nr_vertex + 1, -1); // initializing vector with -1 (assuming no vertex is accessible)
    distances[startV] = 0; // besides the starting one
    queue<int> order_list;
    order_list.push(startV); // starting BFS from the given vertex
    int current;

    while (!order_list.empty()) // there are vertexes to explore
    {
        current = order_list.front();
        for (auto &neighbour: _adjacency[current]) // passing through the vertexes connected to the current one
        {
            if (distances[neighbour] == -1) // if haven't passed through it yet
            {
                distances[neighbour] = distances[current] + 1; //  mark it as accessible
                order_list.push(neighbour);
            }
        }
        order_list.pop(); // finished exploring the current vertex's neighbors
    }

    return distances;
}

void Graph::_dfs(const int& start, const int& marker, vector<int>& mark, stack<int>& exit_time)
{
    /// Solves DFS recursively
    // start = start vertex
    // marker = value to use when finding a visited vertex
    // mark = vector of _nr_vertex + 1 elements; retains the marker of each vertex
    // order_list = vector that determines the longest path in the graph starting with start

    mark[start] = marker; // visited the current vertex
    for (auto &neighbour: _adjacency[start]) // passing through their neighbors
    {
        if (mark[neighbour] == -1) // haven't marked it already
            _dfs(neighbour, marker, mark, exit_time); // continue with the neighbor
    }

    exit_time.push(start);
}

vector<int> Graph::_popEdges(stack<Edge>& st, const Edge& last_edg)
{
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


void Graph::_leveled_dfs(const int& start, vector<int>& parent, vector<int>& level, vector<int>&  return_level, stack<Edge>& expl_edges, vector<vector<int>>& biconex_comps)
{
    if ( parent.size() != _nr_vertex + 1 && level.size() != _nr_vertex + 1 &&  return_level.size() != _nr_vertex + 1 )
    {
        // first iteration
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

        if (parent[child] == -1) // not visited
        {
            expl_edges.push({start, child});
            parent[child] = start;
            level[child] = level[start] + 1;
            return_level[child] = level[child];

            _leveled_dfs(child, parent, level, return_level, expl_edges, biconex_comps);

            return_level[start] = min(return_level[start], return_level[child]);

            if (parent[start] == 0 && nr_children >= 2)
            {
                vector<int> new_comp = _popEdges(expl_edges, {start,child});
                biconex_comps.push_back(new_comp);
            }

            if (parent[start] != 0 && return_level[child] >= level[start])
            {
                vector<int> new_comp = _popEdges(expl_edges, {start,child});
                biconex_comps.push_back(new_comp);
            }
        }

        else if (child != parent[start]) return_level[start] = min(return_level[start], level[child]);
    }
}


/// Procedures for solving the requirements
vector<int> Graph::solveBFS(ifstream &in)
{
    /// Solving BFS from infoarena
    // Warning - the _nr_vertex and _nr_edges must be read before with "in"
    // n = _nr_vertex m = _nr_edges  s = starting node of BFS
    // Out file "bfs.out" - distances from the starting vertex to the others

    int s;
    in >> s;
    _readEdges(in); // Reading the graph's edges
    in.close();
    vector<int> result = _bfs(s);
    return vector<int> (result.begin() + 1, result.end());
}

int Graph::solveDFS(ifstream& in)
{
     /// Solving BFS from infoarena
     /// Using DFS for the non-visited vertexes in order to determine the number of components by
     /// marking them with the same number
    // n = _nr_vertex    m = _nr_edges

    _readEdges(in);
    in.close();
    int result = 0;
    vector<int> components(_nr_vertex + 1, -1);
    stack<int> aux; // not needed in the solving this

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (components[i] == -1) // A new component
        {
            result++;
            _dfs(i, result, components, aux); // mark the other vertex in the component
        }
    }
    return result;
}

vector<int> Graph::solveTopo(ifstream& in)
{
     /// Solving BFS from infoarena
     /// Using DFS for the non-visited vertexes in order to determine the number of components by
     /// marking them with the same number
    // n = _nr_vertex    m = _nr_edges

    _readEdges(in);
    in.close();
    vector<int> components(_nr_vertex + 1, -1);
    vector<int> sol;
    stack<int> aux;

    for (int i = 1; i <= _nr_vertex; ++i)
    {
        if (components[i] == -1) // A new component
        {
            _dfs(i, NULL, components, aux); // mark the other vertex in the component
        }
    }

    while (!aux.empty())
    {
        sol.push_back(aux.top());
        aux.pop();
    }

    return sol;
}

pair<int,vector<vector<int>>> Graph::solveBiconex(ifstream &in)
{
     _readEdges(in);
    in.close();
    vector<vector<int>> sol;
    vector<int> parent;
    vector<int> level;
    vector<int> rtr_lvl;
    stack<Edge> st;

    _leveled_dfs(1, parent, level, rtr_lvl, st, sol);

    vector<int> last_cmp;
    int x,y;
    while (!st.empty())
    {
        x = st.top().from;
        y = st.top().where;
        st.pop();
        last_cmp.push_back(x);
        last_cmp.push_back(y);
    }
    sol.push_back(last_cmp);

    return pair<int, vector<vector<int>>> (sol.size(), sol);
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
        for (int e_idx = 0; e_idx < sol.second[i].size(); ++e_idx)
        {
            out << sol.second[i][e_idx] << " ";
        }
        out << "\n";
    }
}

int main()
{
    infoarenaBiconex();
    return 0;
}
