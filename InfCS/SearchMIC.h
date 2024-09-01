//
// Created by DBL-XQ on 2024/2/27.
//

#ifndef INC_20240227_INFCS_SEARCHMIC_H
#define INC_20240227_INFCS_SEARCHMIC_H


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <ctime>
#include <cmath>

#include <string>
#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <limits>
#include <algorithm>
#include <map>
#include <hash_map>
#include <set>
#include <stack>
#include <queue>
#include <omp.h>
#include <unordered_set>

using namespace std;

class SMIC{

private:

    int graph_size; // the original graph size
    int edge_size; // edge size
    int size_n; // the max node_id+1 in graph

    double d_avg;
    int k_sum;
    int l_sum;
    int k_max;
    int l_max;

public:

    pair<vector<int>, double> MIC;
    std::vector<double> weight;

    std::set<int> nodeset;
    std::unordered_set<int> core_nodeset; // the core_nodeset after some nodes removed

    std::vector<int> k_core;
    std::vector<int> l_core;


    std::vector< std::vector<int> > AdjListOut;
    std::vector< std::vector<int> > AdjListIn;
    // control the storage
    std::vector< __gnu_cxx::hash_map<int, int> > AdjListOut_Loc;
    std::vector< __gnu_cxx::hash_map<int, int> > AdjListIn_Loc;


    std::vector< std::vector<int> > AdjListOut_copy;
    std::vector< std::vector<int> > AdjListIn_copy;

    std::vector< __gnu_cxx::hash_map<int, int> > AdjListOut_Loc_copy;
    std::vector< __gnu_cxx::hash_map<int, int> > AdjListIn_Loc_copy;

    SMIC();
    ~SMIC();

    //// ================================== basic read graph function ============================================
    void readNM(const string attribute_file);
    void readGraph(const string graph_file);
    void readInf(const string inf_file);

    void readsubGraph(const string graph_file, const int p);
    void readsubInf(const string inf_file, const int p);

    void init_loc(const int numV);
    void graphinf();


    void init_Adjlist_copy();
    void init_loc_copy ();
    void clear_copy ();

    string get_dataset(string dataset);
    double avg_inf();

    //// ========================================== basic  function ==================================================
    void add_di_edge(const int start, const int end);
    void remove_di_edge(int start, int end);
    void add_di_copy_edge(const int start, const int end);
    void remove_di_copy_edge(int start, int end);

    vector<pair<int, int>> batch_del_edges(const vector<int> & nodes);
    vector<pair<int, int>> single_del_edges(int node);
    void batch_add_edges(const vector<pair<int, int>> & edges);
    void batch_add_copy_edges(const vector<pair<int, int>> & edges);
    void batch_add_node_copy_edges(const int node, const set<int> curr_set);
    void batch_del_copy_edges(const vector<pair<int, int>> & edges);


    double sum_vector(const vector<int> & vec);
    double avg_vector(const vector<int> & vec);
    double min_vector(const vector<int> & vec);
    double max_vector(const vector<int> & vec);
    int min(int a, int b);
    double find_max_vector(const vector<int> & vec);
    bool containsVector(const vector<std::pair<vector<int>, double>>& L, const vector<int>& X);

    bool compare_adj_sort(const std::pair<int, double>& a, const std::pair<int, double>& b);


    //// ================================== basic compute D-core function ========================================
    vector<int> compute_K0(int graph_size, vector< std::vector<int> > _AdjListIn, vector< std::vector<int> > _AdjListOut);
    vector<int> compute_L0(int graph_size, vector< std::vector<int> > _AdjListIn, vector< std::vector<int> > _AdjListOut);

    void compute_max_kl_core(const int k, const int l);
    vector<pair<int, int>> compute_kl_core(const set<int> nodes, unordered_set<int> &rem_nodes, const int k, const int l, const int s);
    vector<pair<int, int>> compute_copy_kl_core(const set<int> nodes, unordered_set<int> &rem_nodes, const int k, const int l, const int s);
    vector<pair<int, int>> compute_node_kl_core(const set<int> nodes, unordered_set<int> &rem_nodes,  const int max_node, bool & flag, const int k, const int l, const int s);
    vector<pair<int, int>> compute_node_copy_kl_core(const set<int> nodes, unordered_set<int> &rem_nodes,  const int max_node, bool & flag, const int k, const int l, const int s);

    vector<vector<int>> compute_connected_components(unordered_set<int> nodeset, const int k , const int l, const int s);
    vector<vector<int>> compute_connected_components_copy(unordered_set<int> nodeset, const int k , const int l, const int s);
    vector<vector<int>> compute_node_connected_components_copy(unordered_set<int> nodeset, const int max_node, const int k , const int l, const int s);


    void DFS_Scc(int node, set<int> &dfs_rem, vector<pair<int, int>> &res, vector<bool> &visited, const int k, const int l);

    pair<int, double> global_find_min_greedy_node(pair<vector<int>, double> Li, const int k, const int l, const int s);

    pair<vector<int>, double> compute_global_max_comm_v(const pair<vector<int>, double> max_node_L,
                                                 const int max_node, const int k, const int l, const int s);
    pair<vector<int>, double> compute_batch_global_max_comm_v(const pair<vector<int>, double> max_node_L,
                                                        const int max_node, const int k, const int l, const int s);
    pair<vector<int>, double> compute_local_max_comm_v(const pair<vector<int>, double> max_node_L,
                                                        const int max_node, const int k, const int l, const int s);
    pair<vector<int>, double> compute_local_max_comm_v_neigh(const pair<vector<int>, double> max_node_L,
                                                       const int max_node, const int k, const int l, const int s);
    pair<vector<int>, double> compute_local_max_comm(const vector<pair<int, double>>node_sort,
                                                      const int k, const int l, const int s);
    pair<vector<int>, double> compute_local_max_min_comm(const pair<vector<int>, double> max_node_L, const vector<pair<int, double>>node_sort,
                                                     const int k, const int l, const int s);


    //// ======================================= algorithm function ==============================================


    //// -------- basic exact algorithm -------------------

    void basic_random_avg_global_mic(const string graph_file, int k, int l,int s, string seed_mode, string snum);

    void basic_minweight_avg_global_mic(const string graph_file, int k, int l,int s, string seed_mode, string snum);

    void basic_maxweight_avg_local_mic(const string graph_file, int k, int l,int s, string seed_mode, string snum);


    //// ======== improved heuristic algorithm =============
    //global
    void improved_minweight_avg_global_mic(const string graph_file, int k, int l,int s, string seed_mode, string snum, const int p);

    //local
    void improved_maxweight_avg_local_mic(const string graph_file, int k, int l,int s, string seed_mode, string snum, const int p);


    //// ======== min/max/sum global/local algorithm =============
    void improved_minweight_min_global_mic(const string graph_file, int k, int l,int s, string seed_mode, string snum);
    void improved_one_maxweight_min_local_mic(const string graph_file, int k, int l,int s, string seed_mode, string snum);

    void improved_minweight_max_global_mic(const string graph_file, int k, int l,int s, string seed_mode, string snum);
    void improved_maxweight_max_local_mic(const string graph_file, int k, int l,int s, string seed_mode, string snum);

    void improved_sum_mic(const string graph_file, int k, int l,int s, string seed_mode, string snum);



};




#endif //INC_20240227_INFCS_SEARCHMIC_H
