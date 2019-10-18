/**
 *  @file    TemporalDynamics.cpp
 *  @author  Will Oakley
 *  @date    03/06/2018
 *  @version 1.0
 *
 *  @brief Classes to use for studying dynamical processes on temporal random graphs
 *  
 *  @todo split cpp and h
 */

//#include <omp.h>
#include <random>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <tuple>
#include <utility>
#include <numeric>
#include <stack>
#include <fstream>
#include <sstream>

using std::vector;
using std::unordered_set;
using std::unordered_map;
using std::pair;
using std::set;
using std::make_pair;
using std::string;
using std::ofstream;
using std::abs;
using std::stack;
using std::cout;
using std::endl;
using std::stringstream;

/*
 TO DO: 1. CONSIDER OTHER DATA STRUCTURES FOR STORING THE NODE DATA.  HOW TO ALLOW PEOPLE TO
 PLUG IN THEIR OWN DATA STRUCTURE SO ONE CAN DO A DIFFERENT CALCULATION?
 2. MAKE GAMMA NOT HARD CODED IN EXPERIMENT.  HAVE GAMMA FOR INDIVIDUAL NODES AND GAMMA FOR THE WHOLE
 GRAPH.
 */

/**
 *  @brief A class for graph structure
 *  The Graph class contains the adjacency structure of the graph in
 *  the form of an edgelist and auxilliary functions.
 */
class Graph
{
    
public:
    
    /**
     *  @typedef Edge data structure for edges
     *  @typedef Edgelist data structure for edgelist
     *  @typedef EdgeIteratorRange data structure for range of edges
     */
    typedef pair<size_t,size_t> Edge;
    typedef set<Edge> Edgelist;
    typedef pair<Edgelist::iterator,size_t> EdgeIteratorRange;
    
    /**
     *  @return the number of nodes in the graph
     */
    size_t size() const { return sz; }
    
    /**
     *  Compiler generated null constructor
     */
    Graph() = default;
    
    /**
     *  @param el the edgelist
     *  @param sz the number of nodes
     */
    Graph(Edgelist el, size_t sz) : el(el), sz(sz) {}
    
    /**
     *  @return iterator to beginning of the edgelist
     */
    Edgelist::iterator begin()
    {
        return el.begin();
    }
    
    /**
     *  @return iterator to end of the edgelist
     */
    Edgelist::iterator end()
    {
        return el.end();
    }
    
    /**
     *  @return const_iterator to beginning of the edgelist
     */
    Edgelist::const_iterator begin() const
    {
        return el.begin();
    }
    
    /**
     *  @return const_iterator to end of the edgelist
     */
    Edgelist::const_iterator end() const
    {
        return el.end();
    }
    
    /**
     *  @return iterator to beginning of the edgelist
     */
    size_t edges() const {
        return el.size();
    }
    
    /**
     *  @param i a node id
     *  @return in-degree of a node
     */
    size_t in_degree(size_t i) const {
        return range(i).second;
    }
    
    /**
     *  @brief Gets the in-edges of a node.
     *  In-edges of a node i are edges (i,j).  This function returns
     *  an iterator to the first in-edge and the number of in-edges
     *  If a node has no neighbors, a pair containing an iterator
     *  to the first edge after where the node's edges should
     *  have been is returned, along with the value 0.
     *  @param i the node whose in-edges are wanted
     *  @return a pair denoting the in-edge range
     */
    pair<Edgelist::iterator,size_t> range(size_t i) const
    {
        size_t num = 0;
        Edgelist::iterator lb = el.lower_bound(make_pair(i,0));
        while (lb!=end()&&lb->first == i) {
            ++lb;
            ++num;
        }
        return make_pair(el.lower_bound(make_pair(i,0)),num);
    }
    
    /**
     *  Function to erase the first in-edge in an EdgeIteratorRange
     *  and return an EdgeIteratorRange with the next in-edge
     *  and the number of in-edges remaining.
     *  @param p the in-edges of a node
     *  @return the remaining in-edges of a node
     */
    EdgeIteratorRange erase(EdgeIteratorRange p)
    {
        
        return make_pair(el.erase(p.first),--(p.second));
    }
    
    /**
     *  @param e edge to be searched for
     *  @return true if edge is found; false if not
     */
    bool find_edge(const Edge& e) {
        return el.find(e) != el.end();
    }
    
    /**
     *  @param e edge to be inserted
     *  @return iterator to inserted edge or existinge edge and boolean for whether edge was inserted
     */
    pair<Edgelist::iterator,bool> insert_edge(const Edge& e) {
        return el.insert(e);
    }
    
    /**
     * @param e edge to be erased
     * @return the number of edges erased
     */
    size_t erase_edge(const Edge& e) {
        return el.erase(e);
    }
    
private:
    
    Edgelist el; /// @var set of edges
    size_t sz; /// @var sz the number of nodes
};

// really, there should be a "Node Process" class, which defines what kind of process
// we're looking at, and "Updater" would give us certain statistics about that process.

/**
 *  @brief Abstract class to update data held by nodes on graph
 *  This abstract class updates nodes in an inputed Graph
 *  based upon an update process.
 */
class Updater
{
    
public:
    
    /**
     *  @param sz the number of nodes in the Graph
     */
    Updater(size_t sz) : node_data(sz) {}
    
    typedef pair<size_t,size_t> Edge;
    typedef set<Edge> Edgelist;
    typedef pair<Edgelist::iterator,size_t> EdgeIteratorRange;
    typedef stack<EdgeIteratorRange> Edgepath;
    
    /**
     *  @brief a function to update each node's data
     *  Pure virtual function that updates each node's data based upon
     *  some update process.
     *  @param graph the Graph object to use for the update
     *  @return the updated node data
     */
    virtual vector<unordered_map<size_t,size_t>> operator()(Graph&& graph) = 0;
    
    /**
     *  @brief a function to reset the updater
     *  This function is used to reset the node data.
     */
    virtual void reset() = 0;
    
    /**
     *  @brief a function that gets info about the updater
     *  @return a string containing identifying information about the updater
     */
    virtual string info() = 0;
    
protected:
    
    /**
     *  @brief a function to initialize the node data
     *  @param sz the size of the graph
     */
    virtual void node_initializer(size_t sz) = 0;
    
    /**
     *  @brief function to update the node data from two nodes
     *  This function uses the node data from one node to
     *  update the node data from another node.
     *  @param to the node_id who's data to update
     *  @param from the node_id who's data to use for the update
     */
    virtual void update(size_t to, size_t from) = 0;
    
    /**
     *  @brief function that updates the data stored at a node
     *  @param i the node who's node_id to update
     */
    virtual void node_data_update(size_t i) = 0;
    
    Graph graph; /// @var graph on which the update occurs
    vector<unordered_map<size_t,size_t>> node_data; /// @var data structure to store node data
};

class NeighborUpdater : public Updater  {
    
public:
    
    NeighborUpdater(size_t sz, size_t path_len = 1) : Updater(sz), old_node_data(sz), path_len(path_len) {}
    
    /**
     *  @brief a function to update each node's data
     *  Pure virtual function that updates each node's data based upon
     *  some update process.
     *  @param graph the Graph object to use for the update
     *  @return the updated node data
     */
    virtual vector<unordered_map<size_t,size_t>> operator()(Graph&& graph) {
        
        this->graph = graph;
        for (size_t i = 0 ; i < path_len ; ++i) {
            old_node_data = node_data;
            for (Edgelist::iterator it = graph.begin(); it != graph.end() ; it++) {
                update(it->first,it->second);
            }
        }
        
        for (size_t i = 0 ; i < graph.size() ; ++i) {
            node_data_update(i);
        }
        
        return node_data;
    }
    
    /**
     *  @brief a function to reset the updater
     *  This function is used to reset the node data.
     */
    virtual void reset() = 0;
    
    /**
     *  @brief a function that gets info about the updater
     *  @return a string containing identifying information about the updater
     */
    virtual string info() = 0;
    
protected:
    
    /**
     *  @brief a function to initialize the node data
     *  @param sz the size of the graph
     */
    virtual void node_initializer(size_t sz) = 0;
    
    /**
     *  @brief function to update the node data from two nodes
     *  This function uses the node data from one node to
     *  update the node data from another node.
     *  @param to the node_id who's data to update
     *  @param from the node_id who's data to use for the update
     */
    virtual void update(size_t to, size_t from) = 0;
    
    /**
     *  @brief function that updates the data stored at a node
     *  @param i the node who's node_id to update
     */
    virtual void node_data_update(size_t i) = 0;
    
    vector<unordered_map<size_t,size_t>> old_node_data;
    size_t path_len;
};

class FloodingUpdater : public NeighborUpdater  {
    
public:
    
    FloodingUpdater(size_t sz, size_t path_len = 0) : NeighborUpdater(sz,path_len) {
        node_initializer(sz);
    }
    
    /**
     *  @brief a function to reset the updater
     *  This function is used to reset the node data.
     */
    virtual void reset() {
        node_initializer(graph.size());
    }
    
    /**
     *  @brief a function that gets info about the updater
     *  @return a string containing identifying information about the updater
     */
    virtual string info() {
        stringstream ss;
        ss << "Updater=FloodingUpdater;pathLen=" << path_len;
        return ss.str();
    }
    
protected:
    
    /**
     *  @brief a function to initialize the node data
     *  @param sz the size of the graph
     */
    virtual void node_initializer(size_t sz) {
        for (size_t i = 0 ; i < sz ; ++i) {
            node_data[i].clear();
            node_data[i].emplace(i,1);
        }
    }
    
    /**
     *  @brief function to update the node data from two nodes
     *  This function uses the node data from one node to
     *  update the node data from another node.
     *  @param to the node_id who's data to update
     *  @param from the node_id who's data to use for the update
     */
    virtual void update(size_t to, size_t from) {
        node_data[to].insert(old_node_data[from].begin(),old_node_data[from].end());
    }
    
    /**
     *  @brief function that updates the data stored at a node
     *  @param i the node who's node_id to update
     */
    virtual void node_data_update(size_t i) {}
    
};

/**
 *  @brief Abstract class to update node data based on component structure.
 *  This abstract class uses the "reachability" structure of a fixed graph
 *  to define the update procedure.  In other words, the only information
 *  that can be used to update node i are whether paths exist
 *  from other nodes j to node i.  This Updater simultaneously
 *  creates a DAG by modding out by the SCCs in the graph and updates
 *  node data.
 */
class ComponentUpdater : public Updater {
    
public:
    
    /**
     *  @param sz the size of the underlying graph
     */
    ComponentUpdater(size_t sz) : Updater(sz), sccs(sz) {}
    
    /**
     *  @brief function to reset the node data
     *  @see ComponentUpdater::node_initializer(size_t)
     */
    virtual void reset() {
        node_initializer(graph.size());
    }
    
    /**
     *  @brief function to print out info about the updater
     */
    virtual string info() = 0;
    
    /**
     *  @brief function to run update based on a Graph
     *  This call operator updates the graph via a while-loop that
     *  loops while there are still edges available.  It chooses the first
     *  in-edge in the edgelist, and creates an EdgeIteratorRange for the in-edges
     *  of the node who owned the in-edge.  The EdgeIteratorRange is pushed
     *  onto a stack, and a function is called to do a depth first search
     *  style update.  When done, we clean up and redo the algorithm
     *  on the next available edge.
     *  @param graph the graph to update
     *  @return the updated node data
     *  @see ComponentUpdater::advance()
     *  @see ComponentUpdater::scc_update()
     */
    virtual vector<unordered_map<size_t,size_t>> operator()(Graph&& graph)
    {
        this->graph = graph;
        scc_initializer();
        
        while (this->graph.edges() > 0){
            path.push(this->graph.range(find_scc(this->graph.begin()->first)));
            
            visited_sccs.insert(path.top().first->first);
            advance();
            path.pop();
            visited_sccs.clear();
        }
        scc_update();
        
        for (size_t i = 0 ; i < graph.size() ; ++i) {
            node_data_update(i);
        }
        
        return node_data;
    }
    
protected:
    
    /**
     *  @brief this function advances the DFS-style updater
     *  This function attempts to move the DFS-style updater deeper into the graph.
     *  If the updater is at a node without any in-edges, it stops.
     *  If the updater is at a node that hits the scc of a node
     *  belonging to the path, it stops.
     *  If the updater is at a node with available in-edges to follow,
     *  it continues.
     *  @see eval(bool)
     *  @see rewire_edges()
     */
    virtual bool advance()
    {
        if (path.top().second == 0)
        {
            return true;
        }
        if (loop())
        {
            set_scc(path.top().first->first) = find_scc(path.top().first->second);
            path.top() = graph.erase(path.top());
            rewire_edges();
            return false;
        }
        if (path.top().second != 0)
        {
            visited_sccs.insert(find_scc(path.top().first->second));
            path.push(graph.range(find_scc(path.top().first->second)));
        }
        return eval(advance());
    }
    
    /**
     *  @brief this function handles the result of an attempt to advance
     *  If advance was successful, data is taken from the next node.
     *  @see eval(bool)
     *  @see rewire_edges()
     *  @see advance()
     *  @param adv whether the advance was successful
     *  @return whether updater should continue or retreat
     */
    virtual bool eval(bool adv) {
        if (adv)
        {
            path.pop();
            visited_sccs.erase(find_scc(path.top().first->second));
            update(find_scc(path.top().first->first),find_scc(path.top().first->second));
            path.top() = graph.erase(path.top());
            return advance();
        }
        else
        {
            path.pop();
            visited_sccs.erase(path.top().first->second);
            update(find_scc(path.top().first->second),path.top().first->second);
            if (find_scc(path.top().first->first) == find_scc(path.top().first->second))
            {
                size_t node_id = path.top().first->first;
                path.top() = graph.erase(path.top());
                path.top() = graph.range(node_id);
                return advance();
            }
            else
            {
                set_scc(path.top().first->first) = find_scc(path.top().first->second);
                path.top() = graph.erase(path.top());
                rewire_edges();
                return false;
            }
        }
    }
    
    /**
     *  @brief this function checks whether the path has formed a loop
     *  A loop is formed when the path has come back to a node
     *  that was previously in the path.
     *  @return whether there was a loop
     */
    virtual bool loop()
    {
        size_t next = find_scc(path.top().first->second);
        return visited_sccs.find(next) != visited_sccs.end();
    }
    
    /**
     *  @brief this function rewires a node's edges to its scc representative
     *  Each scc has a representative node, and this function takes the edges
     *  of a node in a scc and gives them to the scc representative
     */
    virtual void rewire_edges() {
        while(path.top().second > 0){
            if (find_scc(path.top().first->first) != find_scc(path.top().first->second)) {
                Edge e = make_pair(find_scc(path.top().first->first),find_scc(path.top().first->second));
                graph.insert_edge(e);
            }
            path.top() = graph.erase(path.top());
        }
    }
    
    /**
     *  @brief function to update the node data from two nodes
     *  This function uses the node data from one node to
     *  update the node data from another node.
     *  @param to the node to update
     *  @param from the node to use for the update
     */
    virtual void update(size_t to, size_t from) = 0;
    
    /**
     *  @brief function that updates the data stored at a node
     *  @param i the node who's node_id to update
     */
    virtual void node_data_update(size_t i) = 0;
    
    /**
     * @brief function to initialize the node data
     * @param sz the size of the underlying graph
     */
    virtual void node_initializer(size_t sz) = 0;
    
    /**
     *  @brief function to update the nodes in an scc
     *  This function uses the data held at the representative of an
     *  scc to update all of the other nodes in the scc.
     */
    virtual void scc_update() {
        //        cout << "SCCS" << endl;
        for (size_t i = 0 ; i < graph.size() ; ++i){
            //            cout << i << " " << find_scc(i) << endl;
            update(i,find_scc(i));
        }
    }
    
    /**
     *  @brief function to initialize each node to be it's own scc
     */
    virtual void scc_initializer()
    {
        for (size_t i = 0 ; i < graph.size() ; ++i) {
            sccs[i] = i;
        }
    }
    
    // probably rename this one
    /**
     *  @brief this function allows one to change the scc of a node
     *  @param i the node to query about
     *  @return a reference to where the scc of a node is stored
     */
    virtual size_t& set_scc(size_t i)
    {
        return sccs[i];
    }
    
    /**
     *  @brief this function returns the scc of a node
     *  @param i the node to query about
     *  @return the scc of the node
     */
    virtual size_t find_scc(size_t i) const {
        while (i != sccs[i]) i = sccs[i];
        return i;
    }
    
    vector<size_t> sccs; /// vector storing the sccs of each node
    unordered_set<size_t> visited_sccs; /// hash storing the nodes/sccs the path has visited
    Edgepath path; /// the path taken by the updater
    
};

/**
 *  @brief Class to do a component update based only on reachability
 *  The update in this class only considers whether a node can reach another.
 *  If there is a path from node j to node i, then node i receives
 *  the data stored at node j.  There is no repeated data; each node
 *  has one piece of data and the most information a node can have
 *  about another node is one piece of data.
 */
class ReachabilityUpdater : public ComponentUpdater {
    
public:
    
    /**
     *  @param sz the size of the graph
     *  @see node_intializer(size_t)
     */
    ReachabilityUpdater(size_t sz) : ComponentUpdater(sz) {
        node_initializer(sz);
    }
    
    /**
     *  @brief creates info about the updater
     *  @return a string of information about the updater
     */
    virtual string info() {
        string s = "Updater=Reachability";
        return s;
    }
    
protected:
    
    /**
     *  @brief function defines how data from two nodes is updated
     *  This function simply inserts the data from one node
     *  into another node, without modifying any of the data.
     *  @param to node data to modify
     *  @param from node data to use in the modification
     */
    
    virtual void update(size_t to,size_t from)
    {
        node_data[to].insert(node_data[from].begin(),node_data[from].end());
    }
    
    /**
     *  @brief do nothing
     */
    virtual void node_data_update(size_t i) {}
    
    /**
     *  @brief function that initializes the node data
     *  This function gives each node a single piece of data.
     *  @param sz the size of the graph
     */
    virtual void node_initializer(size_t sz) {
        for (size_t i = 0 ; i < sz ; ++i) {
            node_data[i].clear();
            node_data[i].emplace(i,1);
        }
    }
    
};

/**
 *  @brief This class does a component update using path information
 *  Each node stores the number of distinct temporal paths from each
 *  node on the graph to the node.  A temporal path here does not
 *  take into account the part of the path that lies in each of the
 *  stationary graphs in the temporal graph.  Instead, a temporal
 *  consists of a single node per stationary graph in the temporal
 *  graph, where two consecutive nodes must be connected in the
 *  stationary graph.
 */
class ReachabilityPathUpdater : public ComponentUpdater {
    
public:
    
    /**
     *  @param sz the size of the graph
     *  @see node_initializer(size_t)
     */
    ReachabilityPathUpdater(size_t sz, long long old) : ComponentUpdater(sz), old(old) {
        node_initializer(sz);
    }
    
    /**
     *  @brief gives information about the updater
     *  @return a string giving info about the updater
     */
    virtual string info() {
        stringstream ss;
        ss << "Updater=ReachabilityPath;old=" << old;
        return ss.str();
    }
    
protected:
    
    /**
     *  @brief function defines how data from two nodes is updated
     *  To update node i, this function computes the number of temporal
     *  paths to node i from each reachable node j in the temporal graph
     *  @param to node data to modify
     *  @param from node data to use in the modification
     */
    virtual void update(size_t to,size_t from) {
        //        for (const auto& [key,val] : node_data[from]) {
        //            if (key == to)
        //                node_data[to][key] = 0;
        //            else if (node_data[to].find(key) != node_data[to].end())
        //                node_data[to][key] = node_data[to][key]<val?node_data[to][key]:val;
        //            else
        //                node_data[to][key] = val;
        //        }
        for (auto it = node_data[from].begin() ; it != node_data[from].end() ; it++) {
            if (it->first == to) {
                node_data[to][it->first] = 0;
            } else if (node_data[to].find(it->first) != node_data[to].end()) {
                node_data[to][it->first] = node_data[to][it->first]<it->second?node_data[to][it->first]:it->second;
            } else {
                node_data[to][it->first] = it->second;
            }
        }
    }
    
    /**
     *  @brief do nothing
     */
    virtual void node_data_update(size_t i) {
        for (auto it = node_data[i].begin() ; it != node_data[i].end() ; ) {
            if (it->first != i)
                ++it->second;
            if (old >= 0 && it->second > old) {
                it = node_data[i].erase(it);
            } else {
                ++it;
            }
        }
    }
    
    /**
     *  @brief function that initializes the node data
     *  This function gives each node a single piece of data.
     *  @param sz the size of the graph
     */
    virtual void node_initializer(size_t sz) {
        for (size_t i = 0 ; i < sz ; ++i) {
            node_data[i].clear();
            node_data[i].emplace(i,0);
        }
    }
    
    long long old;
    
};

/**
 *  @brief this abstract class represents an edge-Markovian Temporal Random Graph
 *  An edge-Markovian temporal random graph is defined by two probability
 *  functions on the edges P(e) and Q(e), and an initial distribution.
 */
class eMTRG {
    
public:
    
    typedef pair<size_t,size_t> Edge;
    typedef set<Edge> Edgelist;
    
    /**
     *  @param sz the size of the graph
     */
    eMTRG(size_t sz) : sz(sz), gen(std::random_device()()) {}
    
    /**
     *  This function returns the previous graph and generates
     *  the current graph.  The first time this function is called,
     *  it returns the initial graph, and each further time it is called,
     *  the next graph in the sequence is returned.
     *  @return current element of the sequence.
     *  @see next()
     */
    virtual Graph operator()()
    {
        Graph temp(graph);
        next();
        return temp;
    }
    
    /**
     *  @brief resets the random graph model
     *  This function resets the random graph model to the intial distribution
     *  @see initializer()
     */
    virtual void reset() {
        initializer();
    }
    
    /**
     *  @brief this function returns a string of information
     *  @return a string of info
     */
    virtual string info() = 0;
    
protected:
    
    /**
     *  @brief generates the next random graph
     *  This function generates the next random graph in the
     *  Markov chain.
     *  @see P(Edge)
     *  @see Q(Edge)
     */
    virtual void next()
    {
        //        size_t num_edges = 0;
        for (size_t i = 0 ; i < sz ; ++i)
        {
            for (size_t j = 0 ; j < sz ; ++j) {
                if (i==j) continue;
                Edge e(i,j);
                //                if (graph.find_edge(e))
                //                    ++num_edges;
                std::bernoulli_distribution X(graph.find_edge(e)?1-Q(e):P(e));
                if (X(gen))
                    graph.insert_edge(e);
                else
                    graph.erase_edge(e);
            }
        }
        //        cout << num_edges << endl;
    }
    /**
     *  @brief probability of an edge appearing if it was not present previously
     *  Returns the probabiliy of an edge appearing if that edge
     *  was not present at the previous time step
     *  @param e the edge
     *  @return the probabability
     */
    virtual double P(const Edge& e) = 0;
    
    /**
     *  @brief probability of an edge not appearing if it was present previously
     *  Returns the probabiliy of an edge not appearing if that edge
     *  was present at the previous time step.
     *  @param e the edge
     *  @return the probabability
     */
    virtual double Q(const Edge& e) = 0;
    
    /**
     *  @brief function that generates the intial graph
     */
    virtual void initializer() = 0;
    
    Graph graph;  /// the graph
    size_t sz;  /// size of the graph
    std::mt19937 gen;
};

/**
 *  @brief this class is for generating Erdos-Renyi Temporal Random Graphs
 *  An Erdos Renyi temporal random graph is defined by two parameters
 *  p and alpha.  The functions P(e) and Q(e) are given by
 *  P(e) = alpha*p and Q(e) = alpha(1-p).
 */
class erTRG : public eMTRG {
    
public:
    
    /**
     *  @param sz the size of the graph
     *  @param alpha the correlation
     *  @param p the probability
     *  @see initializer()
     */
    erTRG(size_t sz, double alpha, double p) : eMTRG(sz), p(p), alpha(alpha) {
        initializer();
    }
    
    /**
     *  @brief this function returns a string of information
     *  @return a string of info
     */
    virtual string info() {
        stringstream ss;
        ss << "Graph Model=erTRG;p=" << p << ";alpha=" << alpha;
        return ss.str();
    }
    
protected:
    
    /**
     *  @brief probability of an edge appearing if it was not present previously
     *  Returns the probabiliy of an edge appearing if that edge
     *  was not present at the previous time step.  P(e) = alpha*p
     *  @param e the edge
     *  @return the probability
     */
    virtual double P(const Edge& e) {
        return alpha*p;
    }
    
    /**
     *  @brief probability of an edge not appearing if it was present previously
     *  Returns the probabiliy of an edge not appearing if that edge
     *  was present at the previous time step.  Q(e) = alpha*(1-p)
     *  @param e the edge
     *  @return the probabability
     */
    virtual double Q(const Edge& e) {
        return alpha*(1-p);
    }
    
    /**
     *  @brief function that generates the intial graph
     *  The intial graph is a directed Erdos-Renyi graph
     *  with probability p and no self-edges
     */
    virtual void initializer() {
        Edgelist el;
        std::bernoulli_distribution X(p);
        for (size_t i = 0 ; i < sz ; ++i) {
            for (size_t j = 0 ; j < sz ; ++j) {
                if (i==j) continue;
                if (X(gen))
                    el.insert(make_pair(i,j));
            }
        }
        graph = Graph(el,sz);
    }
    
    double alpha; /// the correlation
    double p; /// the probability
};

//

// TO DO: Should be a combination of a (DYNAMICAL) NODE process and an EDGE process
//           Updater should be in Experiment...

/**
 *  @brief this class represents a Temporal Process, which is a combination
 *  of a Temporal Graph and a Node process
 *  This class contains a Temporal Graph and an Updater.  The Updater
 *  here represents a process by which information is transferred between
 *  nodes.  These two together define a temporal process.
 *  @tparam S the Temporal Random Graph Model
 *  @tparam T the Updater
 */
template <typename S, typename T>
class TemporalProcess {
    
public:
    
    /**
     *  @param emtrg the temporal random graph model
     *  @param updater the updater
     */
    TemporalProcess(S&& emtrg, T&& updater) : emtrg(emtrg), updater(updater) {}
    
    /**
     *  @brief returns the state of the node process
     *  This function applies an updater to a graph and returns
     *  the result of the update
     *  @return the updated node data
     */
    vector<unordered_map<size_t,size_t>> operator()(){
        return updater(emtrg());
    }
    
    /**
     *  @brief resets the TRG and the Updater
     */
    void reset() {
        emtrg.reset();
        updater.reset();
    }
    
    /**
     *  @brief gets info about the temporal process
     *  @return a string of info about the temporal process
     */
    string info() {
        string s = emtrg.info() + "\n" + updater.info();
        return s;
    }
    
protected:
    
    S emtrg; /// the temporal random graph model
    T updater; /// the updater
    
};

/**
 *  @brief this class represents an experiment on a temporal process
 *  This class runs an experiment on a temporal process, where
 *  an experiment is defined by taking certain measurements of
 *  the temporal process.
 */

template <typename S, typename T>
class Experiment {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param filename the output filename for data from measurements
     *  @see measurment(TemporalProcess)
     */
    Experiment(TemporalProcess<S,T> tp, const size_t trials = 10, const string& filename = "data.txt") : tp(tp), trials(trials), fout(filename) {}
    
    void operator()() {
        fout << info() << "\n";
        fout << tp.info() << "\n";
        for (size_t i = 0 ; i < trials ; ++i) {
            cout << "trial " << i+1 << endl;
            // need for loop to run at least twice.
            trial();
            fout << endl;
            tp.reset();
        }
    }
    
protected:
    
    virtual string info() = 0;
    
    virtual void trial() = 0;
    
    /**
     *  @brief the kind of measurement we want to take of the process
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&&) = 0;
    
    ofstream fout; /// filestream for outputting data
    size_t trials; /// the number of disinct trials
    TemporalProcess<S,T> tp;
    
};

/**
 *  @brief this class counts the total amount of data held by the nodes
 *  This class measures the total amount of data held by the nodes
 *  at each iteration of the temporal process.
 */
template<typename S, typename T>
class IterationExperiment : public Experiment<S,T> {
    
public:
   
    /**
     *  @param trials the number of distinct runs of the process
     *  @param iter the amount of discrete time each process runs
     *  @param filename the output filename for data from measurements
     */
    IterationExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const size_t iter = 50, const string& filename = "data.txt") :
    Experiment<S,T>(tp,trials,filename), iter(iter) {}
    
protected:
    
    using Experiment<S,T>::trials;
    using Experiment<S,T>::tp;
    using Experiment<S,T>::fout;
    
    virtual string info() {
        stringstream ss;
        ss << "TRIALS:" << trials << ";Iterations:" << iter << ";Experiment Info:" << specific_info();
        return ss.str();
    }
    
    virtual void trial() {
        for (size_t j = 0 ; j < iter ; ++j) {
            cout << "iteration " << j+1 << endl;
            double m = measurement(tp());
            cout << m << endl;
            fout << m << " ";
        }
    }
    
    virtual string specific_info() = 0;
    
    /**
     *  @brief the kind of measurement we want to take of the process
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&&) = 0;
    
    size_t iter;
};

/**
 *  @brief this class counts the total amount of data held by the nodes
 *  This class measures the total amount of data held by the nodes
 *  at each iteration of the temporal process.
 */
template<typename S, typename T>
class ResidualExperiment : public Experiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param times_identical the number of iterations that result in the same measurement needed to consider the trial over
     *  @param filename the output filename for data from measurements
     */
    ResidualExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const double times_identical = 10, const string& filename = "data.txt") :
    Experiment<S,T>(tp,trials,filename), times_identical(times_identical) {}
    
protected:
    
    using Experiment<S,T>::trials;
    using Experiment<S,T>::tp;
    using Experiment<S,T>::fout;
    
    virtual string info() {
        stringstream ss;
        ss << "TRIALS:" << trials << ";Times Identical:" << times_identical << ";Experiment Info:" << specific_info();
        return ss.str();
    }
    
    virtual void trial() {
        size_t j = 0;
        double m = measurement(tp());
        double mp;
        cout << "iteration " << ++j << endl;
        cout << m << endl;
        fout << m << " ";
        size_t ti = 0;
        do {
            mp = m;
            m = measurement(tp());
            cout << "iteration " << ++j << endl;
            cout << m << endl;
            fout << m << " ";
            if (m == mp) ++ti;
            else ti = 0;
        } while (ti < times_identical);
    }
    
    virtual string specific_info() = 0;
    
    /**
     *  @brief the kind of measurement we want to take of the process
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&&) = 0;
    
    size_t times_identical;
};

/**
 *  @brief this class counts the total amount of data held by the nodes
 *  This class measures the total amount of data held by the nodes
 *  at each iteration of the temporal process.
 */
template<typename S, typename T>
class EqualExperiment : public Experiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param equal_to_value the number of iterations that result in the same measurement needed to consider the trial over
     *  @param filename the output filename for data from measurements
     */
    EqualExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const double equal_to_value = 0, const string& filename = "data.txt") :
    Experiment<S,T>(tp,trials,filename), equal_to_value(equal_to_value) {}
    
protected:
    
    using Experiment<S,T>::trials;
    using Experiment<S,T>::tp;
    using Experiment<S,T>::fout;
    
    virtual string info() {
        stringstream ss;
        ss << "TRIALS:" << trials << ";Equal-to-value:" << equal_to_value << ";Experiment Info:" << specific_info();
        return ss.str();
    }
    
    virtual void trial() {
        size_t j = 0;
        double m = measurement(tp());
        cout << "iteration " << ++j << endl;
        cout << m << endl;
        fout << m << " ";
        while (m != equal_to_value) {
            m = measurement(tp());
            cout << "iteration " << ++j << endl;
            cout << m << endl;
            fout << m << " ";
        }
    }
    
    virtual string specific_info() = 0;
    
    /**
     *  @brief the kind of measurement we want to take of the process
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&&) = 0;
    
    size_t equal_to_value;
};


/**
 *  @brief this class counts the total amount of data held by the nodes
 *  This class measures the total amount of data held by the nodes
 *  at each iteration of the temporal process.
 */
template<typename S, typename T>
class TotalerIterationExperiment : public IterationExperiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param iter the amount of discrete time each process runs
     *  @param filename the output filename for data from measurements
     */
    TotalerIterationExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const size_t iter = 50, const string& filename = "data.txt") :
    IterationExperiment<S,T>(tp,trials,iter,filename) {}
    
protected:
    
    using IterationExperiment<S,T>::iter;
    
    virtual string specific_info() {
        string s = "TotalerIterationExperiment";
        return s;
    }
    
    /**
     *  This function simply gets the total amount of data held by the nodes
     *  @param node_data the hash map of node data
     *  @return the amount of data
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&& node_data)
    {
        double g = 0;
        for (size_t i = 0 ; i < node_data.size() ; ++i) {
            g += node_data[i].size();
        }
        return g;
    }
};

/**
 *  @brief this class counts the total amount of data held by the nodes
 *  This class measures the total amount of data held by the nodes
 *  at each iteration of the temporal process.
 */
template<typename S, typename T>
class TotalerResidualExperiment : public ResidualExperiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param times_identical the number of iterations that result in the same measurement needed to consider the trial over
     *  @param filename the output filename for data from measurements
     */
    TotalerResidualExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const size_t times_identical = 50, const string& filename = "data.txt") :
    ResidualExperiment<S,T>(tp,trials,times_identical,filename) {}
    
protected:
    
    using ResidualExperiment<S,T>::times_identical;
    
    virtual string specific_info() {
        string s = "TotalerResidualExperiment";
        return s;
    }
    
    /**
     *  This function simply gets the total amount of data held by the nodes
     *  @param node_data the hash map of node data
     *  @return the amount of data
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&& node_data)
    {
        double g = 0;
        for (size_t i = 0 ; i < node_data.size() ; ++i) {
            g += node_data[i].size();
        }
        return g;
    }
};

/**
 *  @brief this class counts the total amount of data held by the nodes
 *  This class measures the total amount of data held by the nodes
 *  at each iteration of the temporal process.
 */
template<typename S, typename T>
class TotalerEqualExperiment : public EqualExperiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param equal_to_value the number of iterations that result in the same measurement needed to consider the trial over
     *  @param filename the output filename for data from measurements
     */
    TotalerEqualExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const double equal_to_value = 0, const string& filename = "data.txt") :
    EqualExperiment<S,T>(tp,trials,equal_to_value,filename) {}
    
protected:
    
    using EqualExperiment<S,T>::equal_to_value;
    
    virtual string specific_info() {
        string s = "TotalerEqualExperiment";
        return s;
    }
    
    /**
     *  This function simply gets the total amount of data held by the nodes
     *  @param node_data the hash map of node data
     *  @return the amount of data
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&& node_data)
    {
        double g = 0;
        for (size_t i = 0 ; i < node_data.size() ; ++i) {
            g += node_data[i].size();
        }
        return g;
    }
};

/**
 *  @brief blah blah
 */
template <typename S, typename T>
class ValueTotalerIterationExperiment : public IterationExperiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param iter the amount of discrete time each process runs
     *  @param filename the output filename for data from measurements
     */
    ValueTotalerIterationExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const size_t iter = 50, const string& filename = "data.txt") :
    IterationExperiment<S,T>(tp,trials,iter,filename) {}
    
protected:
    
    using IterationExperiment<S,T>::iter;
    
    virtual string specific_info() {
        string s = "ValueTotalerIterationExperiment";
        return s;
    }
    
    /**
     *  blah
     *  @param node_data the hash map of node data
     *  @return the amount of data
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&& node_data)
    {
        long long total_val = 0;
        for (size_t i = 0 ; i < node_data.size() ; ++i) {
            //            for ([[maybe_unused]] auto [key,val] : node_data[i]){
            //                total_val += val;
            //            }
            for (auto it = node_data[i].begin() ; it != node_data[i].end() ; it++) {
                total_val += it->second;
            }
        }
        return total_val;
    }
};


/**
 *  @brief blah blah
 */
template <typename S, typename T>
class ValueTotalerResidualExperiment : public ResidualExperiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param times_identical the number of iterations that result in the same measurement needed to consider the trial over
     *  @param filename the output filename for data from measurements
     */
    ValueTotalerResidualExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const size_t times_identical = 50, const string& filename = "data.txt") :
    ResidualExperiment<S,T>(tp,trials,times_identical,filename) {}
    
protected:
    
    using ResidualExperiment<S,T>::times_identical;
    
    virtual string specific_info() {
        string s = "ValueTotalerResidualExperiment";
        return s;
    }
    
    /**
     *  blah
     *  @param node_data the hash map of node data
     *  @return the amount of data
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&& node_data)
    {
        long long total_val = 0;
        for (size_t i = 0 ; i < node_data.size() ; ++i) {
            //            for ([[maybe_unused]] auto [key,val] : node_data[i]){
            //                total_val += val;
            //            }
            for (auto it = node_data[i].begin() ; it != node_data[i].end() ; it++) {
                total_val += it->second;
            }
        }
        return total_val;
    }
};

/**
 *  @brief blah blah
 */
template <typename S, typename T>
class ValueTotalerEqualExperiment : public EqualExperiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param equal_to_value the number of iterations that result in the same measurement needed to consider the trial over
     *  @param filename the output filename for data from measurements
     */
    ValueTotalerEqualExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const double equal_to_value = 0, const string& filename = "data.txt") :
    EqualExperiment<S,T>(tp,trials,equal_to_value,filename) {}
    
protected:
    
    using EqualExperiment<S,T>::equal_to_value;
    
    virtual string specific_info() {
        string s = "ValueTotalerEqualExperiment";
        return s;
    }
    
    /**
     *  blah
     *  @param node_data the hash map of node data
     *  @return the amount of data
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&& node_data)
    {
        long long total_val = 0;
        for (size_t i = 0 ; i < node_data.size() ; ++i) {
            //            for ([[maybe_unused]] auto [key,val] : node_data[i]){
            //                total_val += val;
            //            }
            for (auto it = node_data[i].begin() ; it != node_data[i].end() ; it++) {
                total_val += it->second;
            }
        }
        return total_val;
    }
};


/**
 *  @brief blah blah
 */
template <typename S, typename T>
class ValueAveragerIterationExperiment : public IterationExperiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param iter the amount of discrete time each process runs
     *  @param filename the output filename for data from measurements
     */
    ValueAveragerIterationExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const size_t iter = 50, const string& filename = "data.txt") :
    IterationExperiment<S,T>(tp,trials,iter,filename) {}
    
protected:
    
    using IterationExperiment<S,T>::iter;
    
    virtual string specific_info() {
        string s = "ValueAveragerIterationExperiment";
        return s;
    }
    
    /**
     *  blah
     *  @param node_data the hash map of node data
     *  @return the amount of data
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&& node_data)
    {
        double ave = 0;
        size_t num = 0;
        for (size_t i = 0 ; i < node_data.size() ; ++i) {
            //            for ([[maybe_unused]] auto [key,val] : node_data[i]){
            //                ave  = (num*ave + val)/(num+1);
            //                ++num;
            //            }
            for (auto it = node_data[i].begin() ; it != node_data[i].end() ; it++) {
                ave = (num*ave + it->second)/(num+1);
                ++num;
            }
        }
        return ave;
    }
};


/**
 *  @brief blah blah
 */
template <typename S, typename T>
class ValueAveragerResidualExperiment : public ResidualExperiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param times_identical the number of iterations that result in the same measurement needed to consider the trial over
     *  @param filename the output filename for data from measurements
     */
    ValueAveragerResidualExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const size_t times_identical = 50, const string& filename = "data.txt") :
    ResidualExperiment<S,T>(tp,trials,times_identical,filename) {}
    
protected:
    
    using ResidualExperiment<S,T>::times_identical;
    
    virtual string specific_info() {
        string s = "ValueAveragerResidualExperiment";
        return s;
    }
    
    /**
     *  blah
     *  @param node_data the hash map of node data
     *  @return the amount of data
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&& node_data)
    {
        double ave = 0;
        size_t num = 0;
        for (size_t i = 0 ; i < node_data.size() ; ++i) {
            //            for ([[maybe_unused]] auto [key,val] : node_data[i]){
            //                ave  = (num*ave + val)/(num+1);
            //                ++num;
            //            }
            for (auto it = node_data[i].begin() ; it != node_data[i].end() ; it++) {
                ave = (num*ave + it->second)/(num+1);
                ++num;
            }
        }
        return ave;
    }
};

/**
 *  @brief blah blah
 */
template <typename S, typename T>
class ValueAveragerEqualExperiment : public EqualExperiment<S,T> {
    
public:
    
    /**
     *  @param trials the number of distinct runs of the process
     *  @param times_identical the number of iterations that result in the same measurement needed to consider the trial over
     *  @param filename the output filename for data from measurements
     */
    ValueAveragerEqualExperiment(TemporalProcess<S,T> tp, const size_t trials = 10, const size_t equal_to_value = 50, const string& filename = "data.txt") :
    EqualExperiment<S,T>(tp,trials,equal_to_value,filename) {}
    
protected:
    
    using EqualExperiment<S,T>::equal_to_value;
    
    virtual string specific_info() {
        string s = "ValueAveragerEqualExperiment";
        return s;
    }
    
    /**
     *  blah
     *  @param node_data the hash map of node data
     *  @return the amount of data
     */
    virtual double measurement(vector<unordered_map<size_t,size_t>>&& node_data)
    {
        double ave = 0;
        size_t num = 0;
        for (size_t i = 0 ; i < node_data.size() ; ++i) {
            //            for ([[maybe_unused]] auto [key,val] : node_data[i]){
            //                ave  = (num*ave + val)/(num+1);
            //                ++num;
            //            }
            for (auto it = node_data[i].begin() ; it != node_data[i].end() ; it++) {
                ave = (num*ave + it->second)/(num+1);
                ++num;
            }
        }
        return ave;
    }
};

int main(int argc,char* argv[]) {
    
    // parameters for random graph model
    size_t nodes = 1000;
    double alpha = .25;
    double p = .0005;
    // parameters for experiment
    size_t trials = 250;
    size_t iterations = 50;
//    size_t num_identical = 10;
    size_t max_path_len = 10;
//    size_t equal_to_value = 1e6;
    
    // Initialize experiment here.  Using default parameters now.
    for (size_t i = 0 ; i < max_path_len ; ++i) {
        stringstream ss;
        ss << "reachabilityPathExperiment-.0005-old-" << i+1 << ".txt";
        TotalerIterationExperiment<erTRG,ReachabilityPathUpdater> experiment(TemporalProcess<erTRG,ReachabilityPathUpdater>(erTRG(nodes,alpha,p),ReachabilityPathUpdater(nodes,i+1)),trials,iterations,ss.str());
    experiment();
    }
    
    
    // Run experiment
}

