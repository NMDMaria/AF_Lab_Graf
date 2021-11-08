#include <bits/stdc++.h>

using namespace std;

struct Edge
{
    /// Struct that helps with the definition of pairs for the edges of the graphs
    // from = the start vertex of the edge
    // where = the end vertex of the edge
    int from, where;
};

bool pairsortsnddesc(const pair<int,int>& i, const pair<int,int>& j)
{
    /// Function used to sort a vector of pairs descending by the second component
    return i.second > j.second;
}

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
    void _popVertex(const int&, stack<int>&, int&, vector<vector<int>>&, vector<bool>&);
    vector<int> _bfs(const int&);
    void _dfs(const int&, const int&, vector<int>&, stack<int>&);
    void _leveled_dfs(const int&, vector<int>& , vector<int>&, vector<int>& , stack<Edge>&, vector<vector<int>>&, vector<vector<int>>&);
    void _tarjan(const int&, int&, vector<int>&, vector<int>&, stack<int>&, vector<bool>& , int&, vector<vector<int>>&);
    void _HavelHakimi(const int&, const int&, vector<pair<int,int>>& , bool&);
    void _addEdge(const Edge&);
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
    vector<vector<int>> criticalConnections(ifstream&);
    pair<int, vector<vector<int>>> solveCTC(ifstream &);
    void HavelHakimi(ifstream &);
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
    // exit_time = stack in the order of exiting vertex's DFS

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
    // start = the last element that should be removed from stack
    // sol_nr = number of strong connected components
    // sol[i] = vertex's of the strong connected component
    // in[i] = the vertex i was visited in the current component

    sol_nr++; // found a new component
    int aux;
    vector<int> aux2;
    do
    {
        aux = st.top();
        in[aux] = 0; // mark the vertex as not visited anymore
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

    int nr_children = 0; // the start nr. of children in the dfs tree
    for (auto &child: _adjacency[start])
    {
        nr_children ++; // mark as child
        if (parent[child] == -1) // not visited
        {
            expl_edges.push({start, child}); // add the edge we just explored
            parent[child] = start; // mark parent
            level[child] = level[start] + 1;
            return_level[child] = level[child]; // now the lowest level  reached  from the child is his exact level

            _leveled_dfs(child, parent, level, return_level, expl_edges, biconex_comps, critical_edges); // continue DFS tree

            return_level[start] = min(return_level[start], return_level[child]); // passed the rest of the DF tree of the child
            // can modify the lowest level reached in case there was a vertex with a return edge

            if (return_level[child] > level[start]) // the child cant return higher then his parent
            {
                critical_edges.push_back(vector<int>{start,child}); // so this edge is critical in the graph
            }

            if (parent[start] == 0 && nr_children >= 2) // if the root of DF tree has 2 or more  children, means it is a critical vertex
            {
                vector<int> new_comp = _popEdges(expl_edges, {start,child});  // get the biconex component of start
                biconex_comps.push_back(new_comp); // add to solution
            }

            if (parent[start] != 0 && return_level[child] >= level[start]) // the child can return to a lower level then his parent
            {
                // so we reached the end of a biconex component
                vector<int> new_comp = _popEdges(expl_edges, {start,child}); // get th ebiconex component of start
                biconex_comps.push_back(new_comp); // add to solution
            }
        }

        else if (child != parent[start]) return_level[start] = min(return_level[start], level[child]); // update the lowest level reachable with return edges
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
    in_time[start] = time; // mark the discovery time of vertex
    time_return_vert[start] = time; // no return edge known
    in_connection[start] = 1; // mark connected
    connection.push(start); // add to stack
    cout << "Time: " << time << " with vertex " << start;
    cout << "\n In_time: ";
    for (auto &it:in_time) cout << it << " ";
    cout <<"\n Low_in: ";
    for (auto &it:time_return_vert) cout << it << " ";
    cout <<"\n Connected: ";
    for (auto it:in_connection) cout << it << " ";
    cout << "\n";

    for (auto &child: _adjacency[start])
    {
        if (in_time[child] == -1) // not discovered child
        {
            // continue df
            _tarjan(child, time,  in_time, time_return_vert, connection, in_connection, nr_ctc, strong_connected);
            time_return_vert[start] = min(time_return_vert[start], time_return_vert[child]); // update in case a return edge was found in the child's  df sub-tree
        }
        else if(in_connection[child]) time_return_vert[start] = min(time_return_vert[start], in_time[child]);
       // the child was discovered and in the same SCC as parent => return edge exists
       // update in case child was discovered before the last child that updated
    }

    if (time_return_vert[start] == in_time[start]) // a vertex that doesn't have a return edge => end of SCC
    {
         cout <<"\n Time to get a CTC with the starting of DF tree  : " << start << "\n";
        _popVertex(start, connection, nr_ctc, strong_connected, in_connection); // update SCC solution
    }
}

void Graph::_HavelHakimi(const int& n, const int& nr_d, vector<pair<int,int>>& degrees , bool& breakflag)
{
    /// Given a vector of pairs (Vertex, degree) determine if it is a valid graph
    int  m = nr_d/2; // number of edges is sum of degrees/2
    _resize(n, m, 0);
    breakflag = 0;
    vector<vector<bool>> matrix_adj; // access to finding if a edge exists in O(1). but memory + O(n^2)
    matrix_adj.resize(n + 1);
    for (int i = 1; i <= n; ++i) // initialize matrix of adjancency with no edges
        matrix_adj[i].resize(n + 1, 0);
    while (true) //the max nr. of iterations is the number of edges in the supposed graph
    {
        sort(degrees.begin(), degrees.end(), pairsortsnddesc); // Sort the degrees descending
        if (degrees.size() == 0 || degrees[0].second == 0) break; // No vertex left with degree > 0
        pair<int,int> max_dg; // Get the vertex with the highest degree
        max_dg.first = degrees[0].first;
        max_dg.second = degrees[0].second;
        degrees.erase(degrees.begin());
        for (unsigned int k = 0; k < degrees.size() && max_dg.second > 0 ; ++k) // Going through the vertex's
        { // till we have connected the current one with the rest achieving his given degree
            if (!matrix_adj[max_dg.first][degrees[k].first]) // we don't have the current edge in the graph
            {
                matrix_adj[max_dg.first][degrees[k].first] = matrix_adj[degrees[k].first][max_dg.first] = 1; // mark new edge
                degrees[k].second--; // lower degrees
                max_dg.second--;
                if (degrees[k].second < 0) { breakflag = 1; break;} // Connected with a vertex with current degree = 0 => not a possible graph
                m--; // mark the addition of edges

                _addEdge({max_dg.first, degrees[k].first}); // update graph with edge
            }
        }
        if (max_dg.second != 0) {breakflag = 1; break;} // we couldn't find enough vertex's to achieve the degree given => not a graph
    }

    if (m!=0) breakflag = 1; // Couldn't add the number of desired edges so not a valid graph
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
     /// Solving DFS from infoarena
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
     /// Solving Sortare Topologica  from infoarena
     /// Using DF tree and exit times order the nodes so if edge i,j exits i appears before j

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
    // With a leveled DFS find critical edges and find biconex components on an non-oriented conex graph
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

    _leveled_dfs(0, parent, level, rtr_lvl, st, sol, crt_edg); // Graph is conex start from anywhere
    return crt_edg; // crt_edg[i] - crt_edg[i][0] crt_edg[i][1] = critical edge
}

pair <int, vector<vector<int>>> Graph::solveCTC(ifstream &in)
{
    /// Solving CTC from infoarena
    // Using Tarjan's algorithm find the strong connected components of a oriented graph
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
        if (in_time[i] == -1) // for every non visited vertex do Tarjan
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
        in >> aux; // read degree
        if (aux > n) // if degree > number of vertex's => impossible
        {
            breakflag = 1;
            break;
        }
        if (!aux) nr_zeros++;
        m += aux; // the number of edges is the number of degrees/2
        degrees.push_back(make_pair(i, aux));
    }
    if  ( m% 2 == 1 || breakflag ) cout << "Not a graph.\n"; // Sum of edges of non-oriented graph needs to be even. if an odd number of odd degrees -> one left without the number of desired edges
    else   if (nr_zeros == n) cout <<"An empty graph with " << n << " vertex's.\n"; // if all degrees are 0 -> graph with isolated vertex's
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
                    cout << e << " ";
                }
                cout<< endl;
            }
     }
}


/// Functions to solve.
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

int main()
{

    return 0;
}
