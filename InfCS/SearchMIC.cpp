//
// Created by DBL-XQ on 2024/2/27.
//

#include "SearchMIC.h"
#include <stack>
#include <unordered_map>
#include <float.h>
#include <random>

SMIC::SMIC() {
    graph_size = -1; // the original graph size
    edge_size = -1;
    size_n = 0;
    d_avg = 0.0;
    k_sum = 0;
    l_sum = 0;
    k_max = 0;
    l_max = 0;
}

SMIC::~SMIC() {
    weight.clear();
    k_core.clear();
    l_core.clear();
    nodeset.clear();
    core_nodeset.clear();
    for(auto & iter : AdjListOut) {
        iter.clear();
    }
    for(auto & iter : AdjListIn) {
        iter.clear();
    }
    AdjListOut.clear();
    AdjListIn.clear();
    for (auto & iter : AdjListOut_Loc) {
        iter.clear();
    }
    for (auto & iter : AdjListIn_Loc) {
        iter.clear();
    }
    AdjListOut_Loc.clear();
    AdjListIn_Loc.clear();

    for(auto & iter : AdjListOut_copy) {
        iter.clear();
    }
    for(auto & iter : AdjListIn_copy) {
        iter.clear();
    }
    AdjListOut_copy.clear();
    AdjListIn_copy.clear();
    for (auto & iter : AdjListOut_Loc_copy) {
        iter.clear();
    }
    for (auto & iter : AdjListIn_Loc_copy) {
        iter.clear();
    }
    AdjListOut_Loc_copy.clear();
    AdjListIn_Loc_copy.clear();
}

//// ============================= Read Graph and Initial Function ==================================
double SMIC::avg_inf(){
    double sum_inf=0.0;
    for (std::set<int>::iterator it = nodeset.begin(); it != nodeset.end(); ++it) {
        int node = *it;
        sum_inf += weight[node];
    }
    cout << sum_inf << " , "<< nodeset.size()<<endl;
    return sum_inf/size_n;


}


void SMIC::readNM(const string attribute_file){
    cout << "Read attribute ..." << endl;
    std::ifstream infile(attribute_file);
    if (!infile.is_open())
    {
        std::cout << "The file \"" + attribute_file + "\" can NOT be opened\n";
        return;
    }
    infile >> graph_size >> edge_size;
    cout << "graph size: " << graph_size << " edge size: " << edge_size << endl;
    infile.close();
}

void SMIC::readGraph(const string graph_file){
    cout << "Read graph ..." << endl;
    std::ifstream infile(graph_file);
    if (!infile.is_open())
    {
        std::cout << "The file \"" + graph_file + "\" can NOT be opened\n";
        return;
    }
    AdjListOut.resize(graph_size+1);
    AdjListIn.resize(graph_size+1);

    init_loc(graph_size+1);

    for (unsigned int i = 0; i < edge_size; i++){
        unsigned int a, b;
        infile >> a >> b;
        if (a == b)
            continue;
        if (a >= graph_size || b>= graph_size)
            continue;
        //cout << "a: " << a << "\tb: " <<b<< endl;
        add_di_edge(a, b);
        nodeset.insert(a); nodeset.insert(b);
        if(a>=size_n) size_n=a;
        if(b>=size_n) size_n=b;
    }
    size_n ++;
    core_nodeset = {nodeset.begin(), nodeset.end()};
    //copy(nodeset.begin(), nodeset.end(), std::inserter(core_nodeset, core_nodeset.begin()));
    cout << "size_n: "<< size_n<<endl;
    infile.close();

}

void SMIC::readInf(const string inf_file){
    cout << "Read Influenced expectation ..." << endl;
    std::ifstream infile(inf_file);
    if (!infile.is_open())
    {
        std::cout << "The file \"" + inf_file + "\" can NOT be opened\n";
        return;
    }
    for (int i = 0; i < size_n; i++)
        weight.push_back(0.0);

    for (unsigned int i = 0; i < graph_size; i++){
        int node;
        double inf;
        infile >> node >> inf;
        if (node >= graph_size)
            continue;
        //cout << "node: " << node << "\tinf: " <<inf<< endl;
        weight[node] = inf;

    }
    infile.close();
}


void SMIC::readsubGraph(const string graph_file, const int p){
    cout << "Read sub graph ..." << endl;
    std::ifstream infile(graph_file);
    if (!infile.is_open())
    {
        std::cout << "The file \"" + graph_file + "\" can NOT be opened\n";
        return;
    }
    int x = 1/ (p * 0.01);
    AdjListOut.resize((graph_size*x) + 1);
    AdjListIn.resize((graph_size*x) + 1);

    init_loc((graph_size*x) + 1);

    for (unsigned int i = 0; i < edge_size; i++){
        unsigned int a, b;
        infile >> a >> b;
        if (a == b)
            continue;
        if (a >= (graph_size*x) || b>= (graph_size*x))
            continue;
        //cout << "a: " << a << "\tb: " <<b<< endl;
        add_di_edge(a, b);
        nodeset.insert(a); nodeset.insert(b);
        if(a>=size_n) size_n=a;
        if(b>=size_n) size_n=b;
    }
    size_n ++;
    core_nodeset = {nodeset.begin(), nodeset.end()};
    //copy(nodeset.begin(), nodeset.end(), std::inserter(core_nodeset, core_nodeset.begin()));
    cout << "size_n: "<< size_n<<endl;
    infile.close();

}

void SMIC::readsubInf(const string inf_file, const int p){
    cout << "Read Influenced expectation ..." << endl;
    std::ifstream infile(inf_file);
    if (!infile.is_open())
    {
        std::cout << "The file \"" + inf_file + "\" can NOT be opened\n";
        return;
    }
    int x = 1/ (p * 0.01);
    for (int i = 0; i < (graph_size*x) + 1; i++)
        weight.push_back(0.0);

    for (unsigned int i = 0; i < (graph_size*x); i++){
        int node;
        double inf;
        infile >> node >> inf;
        if (node >= (graph_size*x))
            continue;
        //cout << "node: " << node << "\tinf: " <<inf<< endl;
        weight[node] = inf;

    }
    infile.close();
}

//description: This function is used to calculate all nodes's Max K-core-number when there is no out degree constraint;
//Namely, calculate all (k,0)-core for any possible k
vector<int> SMIC::compute_K0(int graph_size, vector< std::vector<int> > _AdjListIn, vector< std::vector<int> > _AdjListOut) {
    vector<int> Kdeg, kBin; //used to calculate all (k,0)-cores
    Kdeg.resize(graph_size);
    kBin.resize(graph_size);
    //Kdeg stores all vertices's max possible K, where there is a (K,0)-core contains it
    //kBin is used to do a bucket sort for all vertices in ascending order of their Kdeg

    vector<int> bin, pos, vert;//used to calculate rows
    bin.resize(graph_size+1);
    pos.resize(graph_size+1);
    vert.resize(graph_size+1);
    //give a certain K,
    //bin is a bucked sort in ascending number of vertices' max L where there is a (K,L)-core contains it
    //pos: which bin a vertex is in
    //vert: the array of all vertices in this row(namely in (K,0)-core)

    int Kmd = -1; //max K degree of the whole graph;
    for (int i = 0; i < graph_size; i++) {
        Kdeg[i] = _AdjListIn[i].size();
        if (Kdeg[i] > Kmd) Kmd = Kdeg[i];
    }
    for (int i = 0; i < graph_size; i++) {
        kBin[Kdeg[i]] += 1;
    }
    int start = 0;       //initialize each bin's start position
    for (int i = 0; i <= Kmd; i++) {
        int num = kBin[i];
        kBin[i] = start;
        start += num;
    }
    //fill bucket's with nodes
    for (int i = 0; i < graph_size; i++) {
        pos[i] = kBin[Kdeg[i]];
        vert[pos[i]] = i;
        kBin[Kdeg[i]] += 1;
    }
    //re-update bin to store each bucket's start position
    for (int i = Kmd; i > 0; i--) {
        kBin[i] = kBin[i - 1];
    }
    kBin[0] = 0;
    // decompose
    for (int i = 0; i < graph_size; i++) {
        int v = vert[i];
        vector<int> vNeighbor = _AdjListOut[v];
        for (int u : vNeighbor) {  //Note there we should use out link to find (v,u) and reduce u's Kdeg
            if (Kdeg[u] > Kdeg[v]) {
                int du = Kdeg[u], pu = pos[u], pw = kBin[du], w = vert[pw];
                pos[u] = pw;    vert[pu] = w;
                pos[w] = pu;    vert[pw] = u;
                kBin[du] += 1;           //pull that higher-order vertex back, and the start position consequently +1
                Kdeg[u] -= 1;
            }
        }
    }
    return Kdeg;
}

vector<int> SMIC::compute_L0(int graph_size, vector< std::vector<int> > _AdjListIn, vector< std::vector<int> > _AdjListOut) {
    vector<int> Ldeg, lBin; //used to calculate all (0,l)-cores
    Ldeg.resize(graph_size);
    lBin.resize(graph_size);

    vector<int> bin, pos, vert;
    bin.resize(graph_size+1);
    pos.resize(graph_size+1);
    vert.resize(graph_size+1);

    int Lmd = -1; //max L degree of the whole graph
    for (int i = 0; i < graph_size; i++) {
        Ldeg[i] = _AdjListOut[i].size();
        if (Ldeg[i] > Lmd) Lmd = Ldeg[i];
    }
    for (int i = 0; i < graph_size; i++) {
        lBin[Ldeg[i]] += 1;
    }
    int start = 0;       //initialize each bin's start position
    for (int i = 0; i <= Lmd; i++) {
        int num = lBin[i];
        lBin[i] = start;
        start += num;
    }
    //fill bucket's with nodes
    for (int i = 0; i < graph_size; i++) {
        pos[i] = lBin[Ldeg[i]];
        vert[pos[i]] = i;
        lBin[Ldeg[i]] += 1;
    }
    //re-update bin to store each bucket's start position
    for (int i = Lmd; i > 0; i--) {
        lBin[i] = lBin[i - 1];
    }
    lBin[0] = 0;
    // decompose
    for (int i = 0; i < graph_size; i++) {
        int v = vert[i];
        vector<int> vNeighbor = _AdjListIn[v];
        for (int u : vNeighbor) {  //Note there we should use out link to find (v,u) and reduce u's Kdeg
            if (Ldeg[u] > Ldeg[v]) {
                int du = Ldeg[u], pu = pos[u], pw = lBin[du], w = vert[pw];
                pos[u] = pw;    vert[pu] = w;
                pos[w] = pu;    vert[pw] = u;
                lBin[du] += 1;           //pull that higher-order vertex back, and the start position consequently +1
                Ldeg[u] -= 1;
            }
        }
    }
    return Ldeg;
}

void SMIC::graphinf() {
    k_core = compute_K0(size_n, AdjListIn, AdjListOut);
    l_core = compute_L0(size_n, AdjListIn, AdjListOut);

    // compute d_avg
    for (int i : nodeset) {
        k_sum += AdjListIn[i].size();
        l_sum += AdjListOut[i].size();
    }
    d_avg = (k_sum + l_sum) * 1.0 / (2.0 * graph_size);
    // compute k_max
    for (unsigned int i = 0; i < k_core.size(); i++) {
        if (k_max < k_core[i]) {
            k_max = k_core[i];
        }
    }
    //compute l_max
    for (unsigned int i = 0; i < l_core.size(); i++) {
        if (l_max < l_core[i]) {
            l_max = l_core[i];
        }
    }
    cout << "d_avg: " << d_avg << " k_max: " << k_max << " l_max: " << l_max<< endl;
}

string SMIC::get_dataset(string dataset){
    int x = 14; // 从第14位开始截取
    string new_str = dataset.substr(x); // 从第 x 位开始截取后面的所有字符
    return new_str;
}
//// ===================================== Basic Function ==========================================

//// ------------------------------------- Edge Function  -------------------------------------
void SMIC::init_loc(const int numV) {
    AdjListOut_Loc.resize(numV);
    AdjListIn_Loc.resize(numV);
    for (int i = 0; i < numV; i++) {
        for (unsigned int j = 0; j < AdjListOut[i].size(); j++) {
            AdjListOut_Loc[i][AdjListOut[i][j]] = int(j);
        }
        for (unsigned int j = 0; j < AdjListIn[i].size(); j++) {
            AdjListIn_Loc[i][AdjListIn[i][j]] = int(j);
        }
    }
}

void SMIC::add_di_edge(const int start, const int end) {
    AdjListOut[start].push_back(end);
    AdjListOut_Loc[start][end] = AdjListOut[start].size() - 1;
    AdjListIn[end].push_back(start);
    AdjListIn_Loc[end][start] = AdjListIn[end].size() - 1;
}

void SMIC::batch_add_edges(const vector<pair<int, int>> & edges){
    for(unsigned int i = 0; i < edges.size(); i++) {
        add_di_edge(edges[i].first, edges[i].second);
        //cout << "add edge: " << edges[i].first<< "-"<< edges[i].second<<endl;
    }
}

void SMIC::add_di_copy_edge(const int start, const int end) {
    AdjListOut_copy[start].push_back(end);
    AdjListOut_Loc_copy[start][end] = AdjListOut_copy[start].size() - 1;
    AdjListIn_copy[end].push_back(start);
    AdjListIn_Loc_copy[end][start] = AdjListIn_copy[end].size() - 1;
}

void SMIC::batch_add_copy_edges(const vector<pair<int, int>> & edges){
    for(unsigned int i = 0; i < edges.size(); i++) {
        add_di_copy_edge(edges[i].first, edges[i].second);
        //cout << "add_copy edge: " << edges[i].first<< "-"<< edges[i].second<<endl;
    }
}

void SMIC::batch_add_node_copy_edges(const int node, const set<int> curr_set){
    vector<pair<int, int>> edges;
    for (int neigh: curr_set){
        auto it_in = find(AdjListIn[node].begin(), AdjListIn[node].end(), neigh);
        if (it_in != AdjListIn[node].end()) {
            edges.push_back(make_pair(neigh, node));
        }

        auto it_out = find(AdjListOut[node].begin(), AdjListOut[node].end(), neigh);
        if (it_out != AdjListOut[node].end()) {
            edges.push_back(make_pair(node, neigh));
        }
    }
    for(unsigned int i = 0; i < edges.size(); i++) {
        add_di_copy_edge(edges[i].first, edges[i].second);
        //cout << "add_copy edge: " << edges[i].first<< "-"<< edges[i].second<<endl;
    }

}


void SMIC::remove_di_edge(const int start, const int end) {
    //cout << "remove edge: " << start << "-"<<end<<endl;
    int idx = AdjListOut_Loc[start][end];
    AdjListOut_Loc[start][AdjListOut[start].back()] = idx;
    swap(AdjListOut[start][idx], AdjListOut[start].back());
    AdjListOut[start].pop_back();

    int idy = AdjListIn_Loc[end][start];
    AdjListIn_Loc[end][AdjListIn[end].back()] = idy;
    swap(AdjListIn[end][idy], AdjListIn[end].back());
    AdjListIn[end].pop_back();
}

void SMIC::remove_di_copy_edge(const int start, const int end) {
    //cout << "remove edge: " << start << "-"<<end<<endl;
    int idx = AdjListOut_Loc_copy[start][end];
    AdjListOut_Loc_copy[start][AdjListOut_copy[start].back()] = idx;
    swap(AdjListOut_copy[start][idx], AdjListOut_copy[start].back());
    AdjListOut_copy[start].pop_back();

    int idy = AdjListIn_Loc_copy[end][start];
    AdjListIn_Loc_copy[end][AdjListIn_copy[end].back()] = idy;
    swap(AdjListIn_copy[end][idy], AdjListIn_copy[end].back());
    AdjListIn_copy[end].pop_back();
}

vector<pair<int, int>> SMIC::batch_del_edges(const vector<int> & nodes){
    int u;
    vector<pair<int, int>> res;
    for(unsigned int i = 0; i < nodes.size(); i++) {
        u = nodes[i];
        while (AdjListOut[u].size() > 0){
            res.push_back(make_pair(u, AdjListOut[u][0]));
            remove_di_edge(u, AdjListOut[u][0]);
        }
        while (AdjListIn[u].size() > 0){
            res.push_back(make_pair(AdjListIn[u][0], u));
            remove_di_edge(AdjListIn[u][0], u);
        }
        AdjListOut[u].clear();
        AdjListIn[u].clear();
    }
    return res;
}

vector<pair<int, int>> SMIC::single_del_edges(const int node) {
    vector<pair<int, int>> res;
    while (AdjListOut[node].size() > 0){
        res.push_back(make_pair(node, AdjListOut[node][0]));
        remove_di_edge(node, AdjListOut[node][0]);
    }
    while (AdjListIn[node].size() > 0){
        res.push_back(make_pair(AdjListIn[node][0], node));
        remove_di_edge(AdjListIn[node][0], node);
    }
    AdjListOut[node].clear();
    AdjListIn[node].clear();

    return res;
}

void SMIC::batch_del_copy_edges(const vector<pair<int, int>> & edges){
    for(unsigned int i = 0; i < edges.size(); i++) {
        remove_di_copy_edge(edges[i].first, edges[i].second);
        //cout << "remove_di_copy_edge: " << edges[i].first<< "-"<< edges[i].second<<endl;
    }
}

//// ------------------------------------- Core Function  -------------------------------------


void SMIC::compute_max_kl_core(const int k, const int l) { //d-order算法
    //cout << "Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    vector<bool> visited (size_n, false);
    queue<int> adjin_q;
    queue<int> adjout_q;
    queue <int> check_q; // 将所有度小于k,l的节点加入队列

    for (auto it = nodeset.begin(); it != nodeset.end(); ++it) {
        int v = *it;
        if (AdjListIn[v].size() < k || AdjListOut[v].size() < l){
            check_q.push(v);
            core_nodeset.erase(v);
        }
        else{
            visited[v] = true;
        }
    }
    double time = 0.0;
    while (!check_q.empty()){
        int node = check_q.front();
        check_q.pop();
        //cout << "node: " << node << "\tcheck_node size: " << check_nodeset.size()<<endl;
        for (auto it = AdjListIn[node].begin(); it != AdjListIn[node].end(); ++it) {
            adjin_q.push(*it);
        }
        while (!adjin_q.empty()){
            int in_neibor = adjin_q.front();
            adjin_q.pop();
            remove_di_edge(in_neibor, node);
            if (visited[in_neibor] == true && (AdjListIn[in_neibor].size() < k || AdjListOut[in_neibor].size() < l)){
                check_q.push(in_neibor); // found in core_nodeset
            }
        }

        for (auto it = AdjListOut[node].begin(); it != AdjListOut[node].end(); ++it) {
            adjout_q.push(*it);
        }
        while (!adjout_q.empty()){
            int out_neibor = adjout_q.front();
            adjout_q.pop();
            remove_di_edge(node, out_neibor);
            if (visited[out_neibor] == true && (AdjListIn[out_neibor].size() < k || AdjListOut[out_neibor].size() < l)){
                check_q.push(out_neibor); // found in core_nodeset
            }
        }

        core_nodeset.erase(node);
        visited[node] = false;
    }
    //cout << "delete time: " << time <<endl;
    cout << "current graph size: " << core_nodeset.size() <<endl;
}


vector<pair<int, int>> SMIC::compute_kl_core(const set<int> nodes, unordered_set<int> &rem_nodes, const int k, const int l, const int s) {
    vector<pair<int, int>> del;
    vector<bool> visited (size_n, false);
    queue <int> check_q;
    queue<int> adjin_q;
    queue<int> adjout_q;

    for (auto it = nodes.begin(); it != nodes.end(); ++it) {
        int v = *it;
        check_q.push(v);
        visited[v] = true;
    }

    while (!check_q.empty()){
        int node = check_q.front();
        check_q.pop();
        //cout << "node: " << node << "\tcheck_node size: " << check_nodeset.size()<< "\tAdjListIn.size: " << AdjListIn[node].size() << "\tAdjListOut.size: " << AdjListOut[node].size()<<endl;
        if (AdjListIn[node].size() < k || AdjListOut[node].size() < l){
            while (AdjListIn[node].size() > 0){
                int in_neibor = AdjListIn[node][0];
                remove_di_edge(in_neibor, node); //AdjListIn has changed!
                del.push_back(make_pair(in_neibor, node));
                 if (visited[in_neibor] == true){
                     check_q.push(in_neibor); // found in core_nodeset
                 }
            }
            while (AdjListOut[node].size() > 0){
                int out_neibor = AdjListOut[node][0];
                remove_di_edge(node, out_neibor);
                del.push_back(make_pair(node, out_neibor));
                if (visited[out_neibor] == true){
                    check_q.push(out_neibor); // found in core_nodeset
                }
            }

            rem_nodes.erase(node);
            visited[node] = false;
        }
    }
    /*for(unsigned int i = 0; i < del.size(); i++) {
        cout <<"del edges: "<< del[i].first<< "-" <<del[i].second<<endl;
    }*/
    return  del;
}


vector<pair<int, int>> SMIC::compute_node_kl_core(const set<int> nodes, unordered_set<int> &rem_nodes, const int max_node, bool &flag, const int k, const int l, const int s) {
    vector<pair<int, int>> del;
    vector<bool> visited (size_n, false);
    queue <int> check_q;

    for (auto it = nodes.begin(); it != nodes.end(); ++it) {
        int v = *it;
        check_q.push(v);
        visited[v] = true;
    }

    while (!check_q.empty()){
        int node = check_q.front();
        check_q.pop();
        if (AdjListIn[node].size() < k || AdjListOut[node].size() < l){
            if (node == max_node){
                flag = false; // (k,l)-core里必须有当前weight最大的节点max_node
                break;
            }
            while (AdjListIn[node].size() > 0){
                int in_neibor = AdjListIn[node][0];
                remove_di_edge(in_neibor, node); //AdjListIn has changed!
                del.push_back(make_pair(in_neibor, node));
                if (visited[in_neibor] == true){
                    check_q.push(in_neibor); // found in core_nodeset
                }
            }
            while (AdjListOut[node].size() > 0){
                int out_neibor = AdjListOut[node][0];
                remove_di_edge(node, out_neibor);
                del.push_back(make_pair(node, out_neibor));
                if (visited[out_neibor] == true){
                    check_q.push(out_neibor); // found in core_nodeset
                }
            }
            rem_nodes.erase(node);
            visited[node] = false;
        }
    }
    /*for(unsigned int i = 0; i < del.size(); i++) {
        cout <<"del edges: "<< del[i].first<< "-" <<del[i].second<<endl;
    }*/
    return  del;
}



vector<pair<int, int>> SMIC::compute_node_copy_kl_core(const set<int> nodes, unordered_set<int> &rem_nodes, const int max_node, bool &flag, const int k, const int l, const int s) {
    vector<pair<int, int>> del;
    vector<bool> visited (size_n, false);
    queue <int> check_q;
    queue<int> adjin_q;
    queue<int> adjout_q;
    int rem_nodes_size = rem_nodes.size();

    double s1 = omp_get_wtime();
    for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
        int v = *it;
        if (AdjListIn_copy[v].size() < k || AdjListOut_copy[v].size() < l){
            if (v == max_node){
                flag = false; // (k,l)-core里必须有当前weight最大的节点max_node
                return del;
            }
            check_q.push(v);
            rem_nodes_size--;
        }
        else{
            visited[v] = true;
        }

        if (rem_nodes_size < s){
            return del;
        }
    }
    //cout << "check_q.size: "<<check_q.size()<<endl;
    while (!check_q.empty()){
        int node = check_q.front();
        check_q.pop();
        visited[node] = false;

        //cout << "node: " << node << "\tcheck_node size: " << check_q.size()<<endl;

        for (auto it = AdjListIn_copy[node].begin(); it != AdjListIn_copy[node].end(); ++it) {
            adjin_q.push(*it);
        }
        while (!adjin_q.empty()){
            int in_neibor = adjin_q.front();
            adjin_q.pop();
            remove_di_copy_edge(in_neibor, node);
            del.push_back(make_pair(in_neibor, node));
            if (visited[in_neibor] == true && (AdjListIn_copy[in_neibor].size() < k || AdjListOut_copy[in_neibor].size() < l)){
                check_q.push(in_neibor); // found in core_nodeset
                if (in_neibor == max_node){
                    flag = false; // (k,l)-core里必须有当前weight最大的节点max_node
                    return del;
                }
            }
        }
        for (auto it = AdjListOut_copy[node].begin(); it != AdjListOut_copy[node].end(); ++it) {
            adjout_q.push(*it);
        }
        while (!adjout_q.empty()){
            int out_neibor = adjout_q.front();
            adjout_q.pop();
            remove_di_copy_edge(node, out_neibor);
            del.push_back(make_pair(node, out_neibor));
            if (visited[out_neibor] == true && (AdjListIn_copy[out_neibor].size() < k || AdjListOut_copy[out_neibor].size() < l)){
                check_q.push(out_neibor); // found in core_nodeset
                if (out_neibor == max_node){
                    flag = false; // (k,l)-core里必须有当前weight最大的节点max_node
                    return del;
                }
            }
        }

        rem_nodes.erase(node);

        if (rem_nodes.size() < s){
            return  del;
        }
    }
//    for(unsigned int i = 0; i < del.size(); i++) {
//        cout <<"del edges: "<< del[i].first<< "-" <<del[i].second<<endl;
//    }
    return  del;
}


vector<pair<int, int>> SMIC::compute_copy_kl_core(const set<int> nodes, unordered_set<int> &rem_nodes,const int k, const int l, const int s) {
    vector<pair<int, int>> del;
    vector<bool> visited (size_n, false);
    queue <int> check_q;
    queue<int> adjin_q;
    queue<int> adjout_q;
    int rem_nodes_size = rem_nodes.size();

    for (auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
        int v = *it;
        if (AdjListIn_copy[v].size() < k || AdjListOut_copy[v].size() < l){
            check_q.push(v);
            rem_nodes_size--;
        }
        else{
            visited[v] = true;
        }

//        if (rem_nodes_size < s){
//            return del;
//        }
    }
    //cout << "check_q.size: "<<check_q.size()<<endl;

    while (!check_q.empty()){
        int node = check_q.front();
        check_q.pop();
        visited[node] = false;

        //cout << "node: " << node << "\tcheck_node size: " << check_q.size()<<endl;

        for (auto it = AdjListIn_copy[node].begin(); it != AdjListIn_copy[node].end(); ++it) {
            adjin_q.push(*it);
        }
        while (!adjin_q.empty()){
            int in_neibor = adjin_q.front();
            adjin_q.pop();
            remove_di_copy_edge(in_neibor, node);
            del.push_back(make_pair(in_neibor, node));
            if (visited[in_neibor] == true && (AdjListIn_copy[in_neibor].size() < k || AdjListOut_copy[in_neibor].size() < l)){
                check_q.push(in_neibor); // found in core_nodeset
            }
        }
        for (auto it = AdjListOut_copy[node].begin(); it != AdjListOut_copy[node].end(); ++it) {
            adjout_q.push(*it);
        }
        while (!adjout_q.empty()){
            int out_neibor = adjout_q.front();
            adjout_q.pop();
            remove_di_copy_edge(node, out_neibor);
            del.push_back(make_pair(node, out_neibor));
            if (visited[out_neibor] == true && (AdjListIn_copy[out_neibor].size() < k || AdjListOut_copy[out_neibor].size() < l)){
                check_q.push(out_neibor); // found in core_nodeset
            }
        }

        rem_nodes.erase(node);

//        if (rem_nodes.size() < s){
//            return  del;
//        }
    }
//    for(unsigned int i = 0; i < del.size(); i++) {
//        cout <<"del edges: "<< del[i].first<< "-" <<del[i].second<<endl;
//    }
    return  del;
}

vector<vector<int>> SMIC::compute_connected_components(unordered_set<int> nodeset, const int k , const int l, const int s) {
    //cout << "Compute strong connected components ..." <<endl;
    vector<int> low_link(size_n, -1);
    vector<int> disc(size_n, -1);
    vector<bool> on_stack(size_n, false);
    stack<int> st;
    int time_count = 0;

    vector<vector<int>> sccResult;

    for (int v : nodeset) {
        if (disc[v] != -1) continue;

        stack<pair<int, vector<int>::iterator>> st_iter;
        st_iter.push({v, AdjListOut[v].begin()});
        disc[v] = time_count;
        low_link[v] = time_count;
        time_count += 1;
        st.push(v);
        on_stack[v] = true;

        while (!st_iter.empty()){
            int node = st_iter.top().first;
            auto& it = st_iter.top().second;

            if (it != AdjListOut[node].end()) {
                int w = *it;
                ++it;
                if (disc[w] == -1) {
                    st_iter.push({w, AdjListOut[w].begin()});
                    disc[w] = time_count;
                    low_link[w] = time_count;
                    time_count += 1;
                    st.push(w);
                    on_stack[w] = true;
                } else if (on_stack[w]) {
                    low_link[node] = min(low_link[node], disc[w]);
                }
            } else {
                // after traversing all the neighbours
                // if it's a root node, pop the stack and print SCC
                if (low_link[node] == disc[node]) {
                    vector<int> scc;
                    while (st.top() != node) {
                        int w = st.top();
                        //cout << w << " ";
                        st.pop();
                        on_stack[w] = false;
                        scc.push_back(w);
                    }
                    //cout << node << "\n";
                    scc.push_back(node);
                    on_stack[node] = false;
                    st.pop();

                    if (scc.size() >= s && scc.size() >= (k+l)){
                        bool flag = true;
                        for(auto node : scc) {
                            if (AdjListIn[node].size() < k || AdjListOut[node].size() < l){
                                //cout << "node: " << node << "\t indeg: " << AdjListIn[node].size()
                                //<< "\t outdeg: " << AdjListOut[node].size() << endl;
                                flag = false; break;
                            }
                        }
                        if (flag == true)
                            sccResult.push_back(scc);
                    }

                }
                st_iter.pop();
                if (!st_iter.empty())
                    low_link[st_iter.top().first] = min(low_link[st_iter.top().first], low_link[node]);
            }
        }
    }

    return sccResult;

}



vector<vector<int>> SMIC::compute_node_connected_components_copy(unordered_set<int> nodeset, const int max_node, const int k , const int l, const int s) {
    //cout << "Compute strong connected components ..." <<endl;
    vector<int> low_link(size_n, -1);
    vector<int> disc(size_n, -1);
    vector<bool> on_stack(size_n, false);
    stack<int> st;
    int time_count = 0;
    bool flag_node;

    vector<vector<int>> sccResult;

    for (int v : nodeset) {
        if (disc[v] != -1) continue;

        stack<pair<int, vector<int>::iterator>> st_iter;
        st_iter.push({v, AdjListOut_copy[v].begin()});
        disc[v] = time_count;
        low_link[v] = time_count;
        time_count += 1;
        st.push(v);
        on_stack[v] = true;

        while (!st_iter.empty()){
            int node = st_iter.top().first;
            auto& it = st_iter.top().second;

            if (it != AdjListOut_copy[node].end()) {
                int w = *it;
                ++it;
                if (disc[w] == -1) {
                    st_iter.push({w, AdjListOut_copy[w].begin()});
                    disc[w] = time_count;
                    low_link[w] = time_count;
                    time_count += 1;
                    st.push(w);
                    on_stack[w] = true;
                } else if (on_stack[w]) {
                    low_link[node] = min(low_link[node], disc[w]);
                }
            }
            else {
                // after traversing all the neighbours
                // if it's a root node, pop the stack and print SCC
                flag_node = false;
                if (low_link[node] == disc[node]) {
                    vector<int> scc;
                    while (st.top() != node) {
                        int x = st.top();
                        st.pop();
                        //cout << x << " ";
                        if (x == max_node){//x == max_node
                            flag_node = true;
                        }
                        on_stack[x] = false;
                        scc.push_back(x);
                    }
                    //cout << node << "\n";
                    if (node == max_node){//x == max_node
                        flag_node = true;
                    }
                    scc.push_back(node);
                    on_stack[node] = false;
                    st.pop();

                    if(flag_node == true){
                        if (scc.size() >= s && scc.size() >= (k+l)){
                            bool flag = true;
                            for(auto node : scc) {
                                if (AdjListIn_copy[node].size() < k || AdjListOut_copy[node].size() < l){
                                    //cout << "node: " << node << "\t indeg: " << AdjListIn_copy[node].size()
                                    //<< "\t outdeg: " << AdjListOut_copy[node].size() << endl;
                                    flag = false; break;
                                }
                            }
                            if (flag == true)
                                sccResult.push_back(scc);
                        }
                    }

                }
                st_iter.pop();
                if (!st_iter.empty())
                    low_link[st_iter.top().first] = min(low_link[st_iter.top().first], low_link[node]);
            }
        }
    }

    return sccResult;
}


vector<vector<int>> SMIC::compute_connected_components_copy(unordered_set<int> nodeset, const int k , const int l, const int s) {
    //cout << "Compute strong connected components ..." <<endl;
    vector<int> low_link(size_n, -1);
    vector<int> disc(size_n, -1);
    vector<bool> on_stack(size_n, false);
    stack<int> st;
    int time_count = 0;

    vector<vector<int>> sccResult;

    for (int v : nodeset) {
        if (disc[v] != -1) continue;

        stack<pair<int, vector<int>::iterator>> st_iter;
        st_iter.push({v, AdjListOut_copy[v].begin()});
        disc[v] = time_count;
        low_link[v] = time_count;
        time_count += 1;
        st.push(v);
        on_stack[v] = true;

        while (!st_iter.empty()){
            int node = st_iter.top().first;
            auto& it = st_iter.top().second;

            if (it != AdjListOut_copy[node].end()) {
                int w = *it;
                ++it;
                if (disc[w] == -1) {
                    st_iter.push({w, AdjListOut_copy[w].begin()});
                    disc[w] = time_count;
                    low_link[w] = time_count;
                    time_count += 1;
                    st.push(w);
                    on_stack[w] = true;
                } else if (on_stack[w]) {
                    low_link[node] = min(low_link[node], disc[w]);
                }
            }
            else {
                // after traversing all the neighbours
                // if it's a root node, pop the stack and print SCC
                if (low_link[node] == disc[node]) {
                    vector<int> scc;
                    while (st.top() != node) {
                        int w = st.top();
                        st.pop();
                        //cout << w << " ";
                        on_stack[w] = false;
                        scc.push_back(w);
                    }
                    //cout << node << "\n";
                    scc.push_back(node);
                    on_stack[node] = false;
                    st.pop();

                    if (scc.size() >= s && scc.size() >= (k+l)){
                        bool flag = true;
                        for(auto node : scc) {
                            if (AdjListIn_copy[node].size() < k || AdjListOut_copy[node].size() < l){
                                //cout << "node: " << node << "\t indeg: " << AdjListIn_copy[node].size()
                                //<< "\t outdeg: " << AdjListOut_copy[node].size() << endl;
                                flag = false; break;
                            }
                        }
                        if (flag == true)
                            sccResult.push_back(scc);
                    }

                }
                st_iter.pop();
                if (!st_iter.empty())
                    low_link[st_iter.top().first] = min(low_link[st_iter.top().first], low_link[node]);
            }
        }
    }

    return sccResult;
}


bool areVectorsEqual(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> sortedA = a;
    std::vector<int> sortedB = b;

    // 对两个 vector 进行排序
    std::sort(sortedA.begin(), sortedA.end());
    std::sort(sortedB.begin(), sortedB.end());

    // 比较排序后的 vector 是否相等
    return sortedA == sortedB;
}

bool SMIC:: containsVector(const std::vector<std::pair<std::vector<int>, double>>& L, const std::vector<int>& X) {
    for (const auto& entry : L) {
        if (areVectorsEqual(entry.first, X)) {
            return true;
        }
    }
    return false;
}


void SMIC::DFS_Scc(int node, set<int> &dfs_rem, vector<pair<int, int>> &res, vector<bool> &visited, const int k, const int l){
    visited[node] = true;
    vector<int> neibor(AdjListOut[node].begin(),AdjListOut[node].end());
    for (int adj : neibor){
        auto it = dfs_rem.find(adj);
        if (it == dfs_rem.end()) { //not in dfs_scc
            continue;
        }
        remove_di_edge(node, adj);
        res.push_back(make_pair(node, adj));
        if (AdjListIn[adj].size() < k || AdjListOut[adj].size() < l){
            if (visited[adj] == false){
                DFS_Scc(adj, dfs_rem, res, visited, k, l);
            }

        }
    }
    //cout << "delete dfs node: " << node << endl;
    dfs_rem.erase(node);
}


double SMIC::sum_vector(const vector<int> & vec) {
    double sum = 0.0;
    for (int i : vec) {
        sum += weight[i];
    }
    return sum;
}


double SMIC::avg_vector(const vector<int> & vec) {
    double sum = 0.0;
    for (int i : vec) {
        sum += weight[i];
    }
    return sum * 1.0 / vec.size();
}

double SMIC::min_vector(const vector<int> & vec) {
    double min = DBL_MAX;
    for (int i : vec) {
        if (weight[i] < min)
            min = weight[i];
    }
    return min;
}

double SMIC::max_vector(const vector<int> & vec) {
    double max = 0.0;
    for (int i : vec) {
        if (weight[i] > max)
            max = weight[i];
    }
    return max;
}


int SMIC::min(int a, int b) {
        return (a < b) ? a : b;
}


bool comparePairs(const pair<vector<int>, double>& a, const pair<vector<int>, double>& b) {
    return a.second > b.second;
}


bool compare_sort(const pair<int, double>& a, const pair<int, double>& b) {
    return a.second > b.second;
}

bool SMIC::compare_adj_sort(const pair<int, double>& a, const pair<int, double>& b) {
    return a.second > b.second;
    if (a.second != a.second) {
        return a.second > a.second; // 按照 double 值降序排列
    }
    else { // 如果 double 值相同，则比较邻节点数量
        int neighcount_a = AdjListOut[a.first].size() + AdjListIn[a.first].size();
        int neighcount_b = AdjListOut[b.first].size() + AdjListIn[b.first].size();
        return neighcount_a > neighcount_b; // 邻节点数量多的排在前面
    }
}


double SMIC::find_max_vector(const vector<int> & vec){
    double max_weight = 0.0;
    for (int i : vec) {
        if (weight[i]>max_weight)
            max_weight = weight[i];
    }
    return max_weight;
}


pair<int, double> SMIC::global_find_min_greedy_node(pair<vector<int>, double> Li, const int k, const int l, const int s){
    pair<int, double> min_greedy_node; min_greedy_node.second = 0.0;
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    vector<pair<int, double>> node_avgweight;
    vector<int> add_vec;
    int delete_node;
    double update_weight;
    for (int i = 0; i< Li.first.size(); i++){
        delete_node = Li.first[i];
        single_del_edges_vec = single_del_edges(delete_node); // delete edge L\ {L[i] \cup delete_node}
        add_vec.assign(Li.first.begin(), Li.first.end());
        set<int>add_set(add_vec.begin(), add_vec.end());
        add_set.erase(delete_node);
        unordered_set<int> rem_set(add_set.begin(), add_set.end());
        single_kl_del_edges_vec = compute_kl_core(add_set,rem_set, k, l, s);
        vector<vector<int>> resc;
        if (rem_set.size() >= s && rem_set.size() >= (k+l)){
            resc = compute_connected_components(rem_set, k, l, s);
        }
        for(const auto& component : resc) {
            update_weight = avg_vector(component);
            if (update_weight > min_greedy_node.second){
                min_greedy_node.first = delete_node;
                min_greedy_node.second = update_weight;
            }
        }
        batch_add_edges(single_del_edges_vec);
        batch_add_edges(single_kl_del_edges_vec);
    }
    return min_greedy_node;

}


pair<vector<int>, double> SMIC::compute_global_max_comm_v(const pair<vector<int>, double> max_node_L, const int max_node, const int k, const int l, const int s){
    pair<vector<int>, double> MIC_v;
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    vector<pair<vector<int>, double>> L_c;
    vector<vector<pair<int, double>>> exist_weight;
    vector<int> add_vec;
    vector<int> d_core_vec;
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    int delete_node;
    double update_weight;
    bool flag_core, flag_scc;
    MIC_v.second = 0.0;

    exist_weight.resize(1);
    for (unsigned int j = 0; j < max_node_L.first.size(); j++) {
        exist_weight[0].push_back(make_pair(max_node_L.first[j], weight[max_node_L.first[j]]));
        d_core_vec.push_back(max_node_L.first[j]);

    }
    L_c.push_back(max_node_L);
    for (unsigned int i = 0; i < L_c.size(); i++) {
        cout << "L_c size is: "<< L_c.size()<< "\tcheck L_c[" << i  <<"] now\tavg_weight: " << L_c[i].second<< endl; //getchar();
        add_vec.assign(L_c[i].first.begin(), L_c[i].first.end());
        set<int> k_core_set(d_core_vec.begin(), d_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(k_core_set.begin(), k_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);

        //sort the nodes in CC[i] in increasing order of weight
        sort(exist_weight[i].begin(), exist_weight[i].end(),
             [](const pair<int, double> & left, const pair<int, double> & right){
                 return left.second < right.second;
             });

        for (unsigned int j = 0; j < L_c[i].first.size(); j++) {
            /*for (int node: max_node_L.first){
                cout << "node: " <<node << "\tindeg: " << AdjListIn[node].size()<< "\toutdeg: " << AdjListOut[node].size()<<endl;
            }*/
            delete_node = exist_weight[i][j].first;
            cout << "delete node: "<<delete_node <<"\t";
            if (delete_node == max_node){
                continue;
            }
            single_del_edges_vec = single_del_edges(delete_node);
            /*if (AdjListIn[max_node].size() < k || AdjListOut [max_node].size() < l){
                batch_add_edges(single_del_edges_vec);
                continue;
            }*/
            set<int>add_set(add_vec.begin(), add_vec.end());
            add_set.erase(delete_node);
            unordered_set<int> rem_set(add_set.begin(), add_set.end());
            flag_core = true;
            single_kl_del_edges_vec = compute_node_kl_core(add_set,rem_set, max_node, flag_core, k, l, s); //找到的kl-core里面必须包含max_node
            vector<vector<int>> resc;
            if (flag_core == true && rem_set.size() >= s && rem_set.size() >= (k+l)){
                resc = compute_connected_components(rem_set, k, l, s);
            }
            for(const auto& component : resc) {
                flag_scc = false;
                if (containsVector(L_c, component) == false){
                    cout << "component.size: "<< component.size() << "\tavg_weight: "<< avg_vector(component) << "\t{ ";
                    for(auto node : component) {
                        if (node == max_node){
                            flag_scc = true; // 新分解出来的scc里面必须包含max_node
                        }
                        //cout << node <<  ' ';
                        cout << node << '-'<< weight[node] << ' ';
                    }cout << "}"<<endl;
                    if (flag_scc == true){
                        update_weight = avg_vector(component);
                        L_c.push_back(make_pair(component, update_weight));
                        vector<pair<int, double>> nodeInfo;
                        for (int node : component){
                            nodeInfo.push_back(make_pair(node, weight[node]));
                        }
                        exist_weight.push_back(nodeInfo);

                        if (update_weight > MIC_v.second){ // keep the max weight
                            MIC_v.first = component;
                            MIC_v.second = update_weight;
                        }
                    }
                }
            }
            batch_add_edges(single_del_edges_vec);
            batch_add_edges(single_kl_del_edges_vec);
        }
        batch_add_edges(batch_del_edges_vec);
        cout <<endl<< "====> current MIC_v size: " << MIC_v.first.size() << "\tinf_avg: " << MIC_v.second << endl;
    }

    return MIC_v;

}


void SMIC::init_Adjlist_copy(){
    AdjListOut_copy.resize(size_n);
    AdjListIn_copy.resize(size_n);
}

void SMIC::init_loc_copy() {
    AdjListOut_Loc_copy.resize(size_n);
    AdjListIn_Loc_copy.resize(size_n);
    for (int i = 0; i < size_n; i++) {
        for (unsigned int j = 0; j < AdjListOut_copy[i].size(); j++) {
            AdjListOut_Loc_copy[i][AdjListOut_copy[i][j]] = int(j);
        }
        for (unsigned int j = 0; j < AdjListIn_copy[i].size(); j++) {
            AdjListIn_Loc_copy[i][AdjListIn_copy[i][j]] = int(j);
        }
    }
}

void SMIC::clear_copy() {
    for(auto & iter : AdjListOut_copy) {
        iter.clear();
    }
    for(auto & iter : AdjListIn_copy) {
        iter.clear();
    }
    //AdjListOut_copy.clear();
    //AdjListIn_copy.clear();
    for (auto & iter : AdjListOut_Loc_copy) {
        iter.clear();
    }
    for (auto & iter : AdjListIn_Loc_copy) {
        iter.clear();
    }
    //AdjListOut_Loc_copy.clear();
    //AdjListIn_Loc_copy.clear();
}

////+++++++++++++++++++++++++++++++++++++++++++++++ local +++++++++++++++++++++++++++++++++++++++++++++++
pair<vector<int>, double> SMIC::compute_local_max_comm_v(const pair<vector<int>, double> max_node_L, const int max_node, const int k, const int l, const int s){
    //cout << "compute_local_max_comm_v..."<<endl;
    //init_Adjlist_copy ();
    //init_loc_copy ();

    pair<vector<int>, double> MIC_v;
    vector<pair<int, int>> single_kl_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    double update_weight;
    bool flag_core, flag_scc;
    MIC_v.second = 0.0;

    queue<int> q;
    int curr_node;
    set<int> cand_set;
    set<int> neigh_set;
    vector<pair<int, double>> neigh_vec;
    vector<pair<int, int>> neighcout_vec;
//    vector<int> Lv_vec = max_node_L.first;
//    set<int> Lv_set (Lv_vec.begin(), Lv_vec.end());
    vector<bool> visited (size_n, false);
    //vector<bool> visited_neigh (size_n, false);
    //vector<int> count_neigh (size_n, 0);
    //vector<int> deg_k (size_n, 0);
    //vector<int> deg_l (size_n, 0);
    //double time1 = 0.0, time2 = 0.0, time3= 0.0;


    int count=0;
    //// iteration
    q.push(max_node);
    visited[max_node] = true;
    while (!q.empty()){
        curr_node = q.front();
        q.pop();
        // add edges (curr_node, cand_set)
        if (q.size() > 0)
            batch_add_node_copy_edges (curr_node, cand_set);

        cand_set.insert(curr_node);

        if (AdjListIn_copy[curr_node].size() >= k && AdjListOut_copy[curr_node].size() >= l && //新加入的点要有用（满足kl条件）
            cand_set.size() >= s && cand_set.size() >= (k+l)){ //AdjListIn_copy[curr_node].size() >= k && AdjListOut_copy[curr_node].size() >= l
            count ++;
            unordered_set<int> rem_set(cand_set.begin(), cand_set.end());
            flag_core = true;

           //double s1 = omp_get_wtime();
            single_kl_del_edges_vec = compute_node_copy_kl_core(cand_set,rem_set, max_node, flag_core, k, l, s); //找到的kl-core里面必须包含max_node
            //double e1 = omp_get_wtime(); time1 += e1-s1;
            //cout << "max_node: "<< max_node<< "\tindeg: "<< AdjListIn_copy[max_node].size() << "\toutdeg: "<< AdjListOut_copy[max_node].size() <<endl;
            //cout << "curr_node: "<< curr_node<< "\tweight: "<< weight[curr_node]<< "\tcand_set.size: "<< cand_set.size()  << "\t rem_set.size: "<< rem_set.size()<< "\tflag: "<< flag_core << "\tcount: "<< count <<endl;

            // prunning rule
            if(true){
                vector<int> rem_vec(rem_set.begin(), rem_set.end());
                if (MIC.second > 0.0){ //not the first MIC_v
                    if (rem_set.size() > MIC.first.size() && (AdjListIn_copy[max_node].size() < k || AdjListOut_copy[max_node].size() < l)){
                        if (avg_vector(rem_vec) < MIC.second){
                            //cout << "Stop it! rem_set.size: "<< rem_set.size() <<"\t weight: " << avg_vector(rem_vec)<<endl;
                            //cout << "time1: "<< time1 << "\ttime2: "<< time2<< "\ttime3: "<< time3 <<endl;
                            return MIC;
                        }
                    }
                }
            }

            vector<vector<int>> resc_vec;
            if (flag_core == true  && rem_set.size() >= s && rem_set.size() >= (k+l)){
                resc_vec = compute_node_connected_components_copy(rem_set, max_node, k, l, s); // 新分解出来的scc里面已经包含max_node
            }
            //double e2 = omp_get_wtime(); time2+=e2-e1;
            if (resc_vec.size() > 0){
                for(const auto& component : resc_vec) {
                    if (MIC_v.second < avg_vector(component)){
                        MIC_v.first = component;
                        MIC_v.second = avg_vector(component);
                    }
                }
                //cout << "Find it! current MIC_v.size: "<< MIC_v.first.size() << "\t avg_weight: "<< MIC_v.second<<endl;
                //cout << "cand_set.size: "<< cand_set.size()  << "\t rem_set.size: "<< rem_set.size() <<endl;
                //cout << "count: "<< count <<endl;
                //cout << "time1: "<< time1 << "\ttime2: "<< time2<< "\ttime3: "<< time3 <<endl;
                return MIC_v;
            }

            batch_add_copy_edges(single_kl_del_edges_vec);
            //double e3 = omp_get_wtime(); time3 += e3-e2;
        }

        //cout << "inserting neighbors"<<endl;
        neigh_set.clear(); ///////////////////////////插入max_node_L中的邻节点
        for (unsigned int i = 0; i < AdjListIn[curr_node].size(); i++){
            int in_neigh = AdjListIn[curr_node][i];
            if (visited[in_neigh] == false){
                neigh_set.insert(in_neigh);
            }
        }

        for (unsigned int i = 0; i < AdjListOut[curr_node].size(); i++){
            int out_neigh = AdjListOut[curr_node][i];
            if (visited[out_neigh] == false){
                neigh_set.insert(out_neigh);
            }
        }

        neigh_vec.clear();
        for (int neigh : neigh_set){
            neigh_vec.push_back(make_pair(neigh, weight[neigh]));
        }
        sort(neigh_vec.begin(), neigh_vec.end(), compare_sort);

        //cout << "curr_node: "<< curr_node<< "\t neigh_vec.size: " << neigh_vec.size() << endl;
        for (unsigned int i = 0; i < neigh_vec.size(); i ++){
            int max_neigh = neigh_vec[i].first;
            if (visited[max_neigh] == false && AdjListIn[max_neigh].size() >= k && AdjListOut[max_neigh].size() >= l){
                //cout << "q.push "<< max_neigh << "\t weight: "<< weight[max_neigh]<<endl; //getchar();
                q.push(max_neigh);
                visited[max_neigh] = true;

            }
        }

    }
    return MIC_v;

}


pair<vector<int>, double> SMIC::compute_local_max_min_comm(const pair<vector<int>, double> max_node_L, const vector<pair<int, double>>node_sort, const int k, const int l, const int s){
    //cout << "compute_local_max_comm_v..."<<endl;
    pair<vector<int>, double> MIC_v; MIC_v.second = 0.0;
    vector<pair<int, int>> single_kl_del_edges_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    double update_weight;

    int curr_node;
    set<int> cand_set;

    //// iteration
    for (int i = 0; i < node_sort.size(); i ++){
        curr_node = node_sort[i].first;
        if (cand_set.size() > 0)
            batch_add_node_copy_edges (curr_node, cand_set);

        cand_set.insert(curr_node);
        if (AdjListIn_copy[curr_node].size() >= k && AdjListOut_copy[curr_node].size() >= l && //新加入的点要有用（满足kl条件）
            cand_set.size() >= s && cand_set.size() >= (k+l)){ //AdjListIn_copy[curr_node].size() >= k && AdjListOut_copy[curr_node].size() >= l
            unordered_set<int> rem_set(cand_set.begin(), cand_set.end());

            single_kl_del_edges_vec = compute_copy_kl_core(cand_set,rem_set, k, l, s);
            cout << "curr_node: "<< curr_node<< "\tweight: "<< weight[curr_node]<< "\tcand_set.size: "<< cand_set.size()  << "\t rem_set.size: "<< rem_set.size() <<endl;

            vector<vector<int>> resc_vec;
            if (rem_set.size() >= s && rem_set.size() >= (k+l)){
                resc_vec = compute_connected_components_copy(rem_set, k, l, s);
            }
            if (resc_vec.size() > 0){
                for(const auto& component : resc_vec) {
                    cout << "component.size: "<< component.size() << " min_inf: "<< min_vector(component)<<endl;
                    if (MIC_v.second <= min_vector(component)){
                        MIC_v.first = component;
                        MIC_v.second = min_vector(component);
                    }
                }
                cout << "Find it! current MIC_v.size: "<< MIC_v.first.size() << "\t min_weight: "<< MIC_v.second<<endl;
                cout << "cand_set.size: "<< cand_set.size()  << "\t rem_set.size: "<< rem_set.size() <<endl;
                //cout << "count: "<< count <<endl;
                return MIC_v;
            }

            batch_add_copy_edges(single_kl_del_edges_vec);
        }

    }
    return MIC_v;

}


//// ================================ Search MIC Algorithm =============================================

//// ======== basic exact algorithm =============

void SMIC::basic_random_avg_global_mic(const string dataset, int k, int l,int s, string seed_mode, string snum){
    cout << "====================== Naive random avg MIC begin ... ======================" << endl;
    double startTime, endTime;
    double coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> L; //int: connetcted component double: avg weight
    vector<pair<vector<int>, double>> output_res;
//    map<int, int> add_vec_loc;
//    map<int, int> del_nodes_order_vec_loc;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> d_core_vec; // store all the nodes in the sccs
    vector<int> del_nodes_vec;
    vector<int> del_nodes_order_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_iter_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    vector<double> weight_vec;
    double update_weight; // updated weight
    int delete_node;
    double avg_vec;


    startTime = omp_get_wtime();
    coreStartTime = omp_get_wtime();
    cout << "=====> Step1: Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    compute_max_kl_core(k,l); // compute the maximum (k,l)-core in G
    coreEndTime = omp_get_wtime();
    cout << "=====> Step2: Compute strong connected components ..." <<endl;
    res = compute_connected_components(core_nodeset, k, l, s); // compute all strong connected components (scc) on the maximum (k,l)-core

    //// construct initial L0
    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]); // compute the average weight of each sCC
        if (res[i].size() < (k+l) || res[i].size() < s) // make sure nodes in each sCC at least larger than k = each node has at least k neighbors
            continue;

        L.push_back(make_pair(res[i], avg_vec));
        output_res.push_back(make_pair(res[i], avg_vec));
        cout << "res size: " << res[i].size() << "\tres avg: " << avg_vec << endl;
    }
    //// store all the nodes in the maximum d-core (d-core)
    for (unsigned int i = 0; i < L.size(); i++) {
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            d_core_vec.push_back(L[i].first[j]);
        }
    }
    cout << "Init number of connect components of maximal D-core: " << output_res.size() << endl;

    cout << "=====> Step3: Decomposition begin ..." << endl;
    for (unsigned int i = 0; i < L.size(); i++) {
        if (L[i].first.size() < s || L[i].first.size()  < (k+l)) continue;
        add_vec.assign(L[i].first.begin(), L[i].first.end());
        del_nodes_order_vec.assign(L[i].first.begin(), L[i].first.end());
        random_shuffle(del_nodes_order_vec.begin(), del_nodes_order_vec.end()); //// random

//        /// store the location of the delete element
//        del_nodes_order_vec_loc.clear();
//        for (unsigned int j = 0; j < del_nodes_order_vec.size(); j++) {
//            del_nodes_order_vec_loc[del_nodes_order_vec[j]] = j;
//        }
//        /// store the location of element in L[i]
//        add_vec_loc.clear();
//        for (unsigned int j = 0; j < add_vec.size(); j++) {
//            add_vec_loc[add_vec[j]] = j;
//        }
        set<int> d_core_set(d_core_vec.begin(), d_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end()); // add_set = nodes in L[i]
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(d_core_set.begin(), d_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec); // delete edge L\L[i]

        //// delete node in L[i]
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            /*for (int node: nodeset){
                cout << "node: " <<node << "\tindeg: " << AdjListIn[node].size()<< "\toutdeg: " << AdjListOut[node].size()<<endl;
            }*/
            set<int>add_set(add_vec.begin(), add_vec.end());
            delete_node = del_nodes_order_vec[j];

            cout << "delete node: " << delete_node << "\t";

            single_del_edges_vec = single_del_edges(delete_node); // delete edge L\ {L[i] \cup delete_node}

            add_set.erase(delete_node);
            unordered_set<int> rem_set(add_set.begin(), add_set.end());
            single_kl_del_edges_vec = compute_kl_core(add_set,rem_set, k, l, s);
            vector<vector<int>> resc;
            if (rem_set.size() >= s && rem_set.size() >= (k+l)){
                resc = compute_connected_components(rem_set, k, l, s);
            }
            for(const auto& component : resc) {
                cout << "component.size: "<< component.size() << "\tavg_weight: "<< avg_vector(component) << "\t{ ";
                for(auto node : component) {
                    cout << node << ' ';
                }cout << "}"<<endl;
                if (containsVector(L, component) == false){
                    update_weight = avg_vector(component);
                    L.push_back(make_pair(component, update_weight));
                    output_res.push_back(make_pair(component, update_weight));
                }
            }
            batch_add_edges(single_del_edges_vec);
            batch_add_edges(single_kl_del_edges_vec);
        }
        batch_add_edges(batch_del_edges_vec);
    }

    cout << "naive random avg MIC end" << endl;

    //// sort the CC in output_res in descending order of avg_weight
    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    for (unsigned int i = 0; i < output_res.size(); i++) {
        cout << i << " inf_avg is: " << output_res[i].second << " size is: " << output_res[i].first.size() << "\t{ ";
            for(auto node : output_res[i].first) {
                cout << node << ' ';
            }
            cout << "}"<<endl;
    }

    // output result
    MIC = output_res[0];
    cout << "=====> Step4: MIC avg_inf is: "<< MIC.second << " size: " << MIC.first.size() << "\t{ ";
    for (int i = 0; i < MIC.first.size(); i ++){
        cout << MIC.first[i] << " ";
    }
    cout << "}"<< endl;

    endTime = omp_get_wtime();
    cout << "k: " << k << " l: " << l << " s: " << s << endl;
    cout << "startTime: " << startTime << "\t endTime: " << endTime<< endl;
    cout << "Total Time : " << (double)(endTime - startTime)  << " s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;

    // write result to file
    string k_name = "k_" + std::to_string(k) + "_";
    string l_name = "l_" + std::to_string(l) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    string filename = dataset + "/result/naive_random_avg_topr_" + k_name + l_name + s_name + ".txt";
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime)  << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime)  << endl;
    outfile << MIC.second << " " << MIC.first.size() <<"\t";
    for (int i = 0; i < MIC.first.size(); i ++){
        outfile << MIC.first[i] << " ";
    }
    outfile << endl;

    outfile.close();

}

void SMIC::basic_minweight_avg_global_mic(const string dataset, int k, int l,int s, string seed_mode, string snum){
    cout << "====================== Naive greedy avg MIC begin ... ======================" << endl;
    double startTime, endTime;
    double coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> L; //int: connetcted component double: avg weight
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> d_core_vec; // store all the nodes in the sccs
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_iter_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    vector<double> weight_vec;
    vector<vector<pair<int, double>>> exist_weight;
    double update_weight; // updated weight
    int delete_node;
    double avg_vec;

    startTime = omp_get_wtime();

    cout << "=====> Step1: Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    coreStartTime = omp_get_wtime();
    compute_max_kl_core(k,l); // compute the maximum (k,l)-core in G
    coreEndTime = omp_get_wtime();

    cout << "=====> Step2: Compute strong connected components ..." <<endl;
    res = compute_connected_components(core_nodeset, k, l, s); // compute all strong connected components (scc) on the maximum (k,l)-core

    //// construct L0
    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]);
        if (res[i].size() < (k+l) || res[i].size() < s)
            continue;

        L.push_back(make_pair(res[i], avg_vec));
        output_res.push_back(make_pair(res[i], avg_vec));
        cout << "res size: " << res[i].size() << "\tres avg: " << avg_vec << endl;
    }

    //// store all the nodes in the maximum d-core (d-core)
    exist_weight.resize(L.size());
    for (unsigned int i = 0; i < L.size(); i++) {
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            exist_weight[i].push_back(make_pair(L[i].first[j], weight[L[i].first[j]]));
            d_core_vec.push_back(L[i].first[j]);
        }
    }
    cout << "Init number of connect components of maximal D-core: " << output_res.size() << endl;

    cout << "=====> Step3: Decomposition begin ..." << endl;
    for (unsigned int i = 0; i < L.size(); i++) {
        if (L[i].first.size() < s || L[i].first.size()  < (k+l))
            continue;
        add_vec.assign(L[i].first.begin(), L[i].first.end());
        set<int> d_core_set(d_core_vec.begin(), d_core_vec.end());
        set<int> add_set(add_vec.begin(), add_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(d_core_set.begin(), d_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec);

        //// greedy
        sort(exist_weight[i].begin(), exist_weight[i].end(),
             [](const pair<int, double> & left, const pair<int, double> & right){
                 return left.second < right.second;
             });

        batch_iter_del_edges_vec.clear();
        //// delete node in L[i]
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            /*for (int node: nodeset){
                cout << "node: " <<node << "\tindeg: " << AdjListIn[node].size()<< "\toutdeg: " << AdjListOut[node].size()<<endl;
            }*/

            delete_node = exist_weight[i][j].first;
            cout << "delete node: " << delete_node << "\t";
            single_del_edges_vec = single_del_edges(delete_node); // delete edge L\ {L[i] \cup delete_node}

            set<int>add_set(add_vec.begin(), add_vec.end());
            add_set.erase(delete_node);
            unordered_set<int> rem_set(add_set.begin(), add_set.end());
            single_kl_del_edges_vec = compute_kl_core(add_set,rem_set, k, l, s);
            vector<vector<int>> resc;
            if (rem_set.size() >= s && rem_set.size() >= (k+l)){
                resc = compute_connected_components(rem_set, k, l, s);
            }
            for(const auto& component : resc) {
                cout << "component.size: "<< component.size() << "\tavg_weight: "<< avg_vector(component) << "\t{ ";
                for(auto node : component) {
                    cout << node << ' ';
                }cout << "}"<<endl;
                if (containsVector(L, component) == false){
                    update_weight = avg_vector(component);
                    L.push_back(make_pair(component, update_weight));

                    vector<pair<int, double>> nodeInfo;
                    for (int node : component){
                        nodeInfo.push_back(make_pair(node, weight[node]));
                    }
                    exist_weight.push_back(nodeInfo);

                    output_res.push_back(make_pair(component, update_weight));
                }
            }
            batch_add_edges(single_del_edges_vec);
            batch_add_edges(single_kl_del_edges_vec);
        }
        batch_add_edges(batch_del_edges_vec);
    }

    cout <<endl << "Naive greedy avg MIC end" << endl;

    //// sort the CC in output_res in descending order of avg_weight
    sort(output_res.begin(), output_res.end(),
         [](const pair<vector<int>, double> & left, const pair<vector<int>, double> & right){
             return left.second > right.second;
         });

    for (unsigned int i = 0; i < output_res.size(); i++) {
        cout << i << " inf_avg is: " << output_res[i].second << " size is: " << output_res[i].first.size() << "\t{ ";
        for(auto node : output_res[i].first) {
            cout << node << ' ';
        }
        cout << "}"<<endl;
    }

    // output result
    MIC = output_res[0];
    cout << "=====> Step4: MIC avg_inf is: "<< MIC.second << " size: " << MIC.first.size() << "\t{ ";
    for (int i = 0; i < MIC.first.size(); i ++){
        cout << MIC.first[i] << " ";
    }
    cout << "}"<< endl;

    endTime = omp_get_wtime();
    cout << "k: " << k << " l: " << l << " s: " << s << endl;
    cout << "startTime: " << startTime << "\t endTime: " << endTime<< endl;
    cout << "Total Time : " << (double)(endTime - startTime)  << " s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;

    // write result to file
    string k_name = "k_" + std::to_string(k) + "_";
    string l_name = "l_" + std::to_string(l) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    string filename = dataset + "/result/naive_greedy_avg_topr_" + k_name + l_name + s_name + ".txt";
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime)  << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime)  << endl;
    outfile << MIC.second << " " << MIC.first.size() <<"\t";
    for (int i = 0; i < MIC.first.size(); i ++){
        outfile << MIC.first[i] << " ";
    }
    outfile << endl;

    outfile.close();

}

void SMIC::basic_maxweight_avg_local_mic(const string dataset, int k, int l,int s, string seed_mode, string snum){
    cout << "====================== Improved greedy avg local MIC begin ... ======================" << endl;
    double startTime, endTime;
    double coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> L;
    vector<int> add_vec; // after delete one node, new vector
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    double update_weight; // updated weight
    int delete_node;
    double avg_vec;
    double p =0.0;//p =0.005


    vector<pair<int, double>> nodeset_sort;
    vector<int> node_L_loc (size_n, -1);
    pair<vector<int>, double> MIC_v;


    startTime = omp_get_wtime();
    cout << "=====> Step1: Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    coreStartTime = omp_get_wtime();
    compute_max_kl_core(k,l); // compute the maximum (k,l)-core in G
    coreEndTime = omp_get_wtime();
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;
    getchar();

    cout << "=====> Step2: Compute strong connected components ..." <<endl;
    res = compute_connected_components(core_nodeset, k, l, s); // compute all strong connected components (scc) on the maximum (k,l)-core

    //// construct initial L0
    MIC.second = 0.0;
    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]); // compute the average weight of each sCC
        if (res[i].size() < (k+l) || res[i].size() < s) // make sure nodes in each sCC at least larger than k = each node has at least k neighbors
            continue;

        L.push_back(make_pair(res[i], avg_vec));
        for (int j = 0; j < res[i].size(); j ++){
            int node = res[i][j];
            node_L_loc[node] = L.size() - 1;
        }

        cout << "res size: " << res[i].size() << "\tres avg: " << avg_vec << endl;

        if (avg_vec > MIC.second){
            MIC.first = res[i];
            MIC.second = avg_vec;
        }
    }
    cout << "=====> MIC size: " << MIC.first.size() << "\tinf_avg: " << MIC.second << endl;
    cout << "Init number of connect components of maximal D-core: " << L.size() << endl;

    cout << "=====> Step3: Sort node in V in decreasing order of weight ..." << endl;
    for (int node: nodeset){
        nodeset_sort.push_back(make_pair(node,weight[node]));
    }
    sort(nodeset_sort.begin(), nodeset_sort.end(), compare_sort);


    cout << "=====> Step4: Decomposition begin ..." << endl;
    for (int i = 0; i < nodeset_sort.size(); i ++){
        int max_node = nodeset_sort[i].first;
        cout << "check node: " << max_node << "\t weight: "<< weight[max_node] << "\t L id is: "<< node_L_loc[max_node]<< endl; //getchar();
        if (weight[max_node] <= MIC.second){
            break; //若v的weight小于当前最大avg，由于在所有包含v的社区中，v的weight是最大的，因此合并其他weight更小的节点只会使avg减小
        }
        int max_node_Lid = node_L_loc[max_node];
        if (max_node_Lid == -1) // max_node_Lid == -1 <==> AdjListIn[max_node] < k || AdjListOut[max_node] < l
            continue;
        pair<vector<int>, double> max_node_L;
        max_node_L.first = L[max_node_Lid].first;
        max_node_L.second = L[max_node_Lid].second;
        if (max_node_L.first.size() < s || max_node_L.first.size()  < (k+l))
            continue;

        ////////////////// 找到包含max_node的所有社区中 avg influenced expectation 最大的社区 //////////////////
        cout << "=====> Step4.1: Find the maximum weight comm including node "<< max_node << endl;
        MIC_v = compute_global_max_comm_v(max_node_L, max_node, k, l, s);
        cout << "=====> Step4.2: MIC_v size: "<< MIC_v.first.size() << "\tinf_avg: " << MIC_v.second << endl;
        if (MIC_v.second > MIC.second){
            MIC.first = MIC_v.first;
            MIC.second = MIC_v.second;
        }
        cout << "=====> MIC size: " << MIC.first.size() << "\tinf_avg: " << MIC.second << endl;


        // 删除当前max_node，分解原本的L[max_node_Lid]，生成新的D-core scc, 并更新其余节点的max_node_Lid
        cout << "=====> Step4.2: Delete current max_node "<< max_node << endl;
        delete_node = max_node;
        node_L_loc[max_node] = -1;
        //cout << "delete node: " << delete_node << "\t";
        single_del_edges_vec = single_del_edges(delete_node); // delete edge L\ {L[i] \cup delete_node}

        add_vec.assign(max_node_L.first.begin(), max_node_L.first.end());
        set<int>add_set(add_vec.begin(), add_vec.end());
        add_set.erase(delete_node);
        unordered_set<int> rem_set(add_set.begin(), add_set.end());
        single_kl_del_edges_vec = compute_kl_core(add_set,rem_set, k, l, s);
        vector<vector<int>> resc;
        if (rem_set.size() >= s && rem_set.size() >= (k+l)){
            resc = compute_connected_components(rem_set, k, l, s);
        }
        set<int> ori_set(add_vec.begin(), add_vec.end());
        set<int> fin_set;
        set<int> del_set;

        for(const auto& component : resc) {
            if (containsVector(L, component) == false){
                // cout << "component.size: "<< component.size() << "\tavg_weight: "<< avg_vector(component) << "\t{ ";
                /*for(auto node : component) {
                    cout << node <<  ' ';
                    // cout << node << '-'<< weight[node] << ' ';
                }cout << "}"<<endl;*/
                update_weight = avg_vector(component);
                L.push_back(make_pair(component,update_weight));
                for (int node: component) {
                    node_L_loc[node] = L.size() - 1; //update location L of remain nodes
                    fin_set.insert(node);
                }
            }
        }
        set_difference(ori_set.begin(), ori_set.end(),
                       fin_set.begin(), fin_set.end(),
                       inserter(del_set, del_set.begin()));
        for (int node : del_set){
            node_L_loc[node] = -1;
        }

    }

    cout <<endl << "Improved greedy avg local MIC end" << endl;

    // output result
    cout << "=====> Step5: MIC avg_inf is: "<< MIC.second << " size: " << MIC.first.size() << "\t{ ";
    for (int i = 0; i < MIC.first.size(); i ++){
        cout << MIC.first[i] << " ";
    }
    cout << "}"<< endl;

    endTime = omp_get_wtime();
    cout << "k: " << k << " l: " << l << " s: " << s << endl;
    cout << "startTime: " << startTime << "\t endTime: " << endTime<< endl;
    cout << "Total Time : " << (double)(endTime - startTime)  << " s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;

    // write result to file
    string k_name = "k_" + std::to_string(k) + "_";
    string l_name = "l_" + std::to_string(l) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    string filename = dataset + "/result/improved_greedy_avg_local_mic_" + k_name + l_name + s_name + ".txt";
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime)  << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime)  << endl;
    outfile << MIC.second << " " << MIC.first.size() <<"\t";
    for (int i = 0; i < MIC.first.size(); i ++){
        outfile << MIC.first[i] << " ";
    }
    outfile << endl;

    outfile.close();

}




//// ======== improved heuristic algorithm =============
////+++++++++++++++++++++++++++++++++++++++++++++++ global +++++++++++++++++++++++++++++++++++++++++++++++
void SMIC::improved_minweight_avg_global_mic(const string dataset, int k, int l,int s, string seed_mode, string snum, const int p){
    cout << "====================== Improved greedy avg global MIC begin ... ======================" << endl;
    double startTime, endTime;
    double coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> L; //int: connetcted component double: avg weight
    vector<pair<vector<int>, double>> res_L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> d_core_vec; // store all the nodes in the sccs
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_iter_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    vector<double> weight_vec;
    vector<vector<pair<int, double>>> exist_weight;
    double update_weight; // updated weight
    int delete_node;
    double avg_vec;

    startTime = omp_get_wtime();

    cout << "=====> Step1: Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    coreStartTime = omp_get_wtime();
    compute_max_kl_core(k,l); // compute the maximum (k,l)-core in G
    coreEndTime = omp_get_wtime();

//    double sum_inf = 0.0;
//    double avg_inf;
//    for (std::unordered_set<int>::iterator it = core_nodeset.begin(); it != core_nodeset.end(); ++it) {
//        int node = *it;
//        sum_inf += weight[node];
//    }
//    avg_inf=sum_inf/core_nodeset.size();
//    cout << "sum inf: "<< sum_inf << "avg inf: "<< avg_inf<<endl;

    cout << "=====> Step2: Compute strong connected components ..." <<endl;
    res = compute_connected_components(core_nodeset, k, l, s); // compute all strong connected components (scc) on the maximum (k,l)-core

    //// construct initial L0
    MIC.second = 0.0;
    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]); // compute the average weight of each sCC
        if (res[i].size() < (k+l) || res[i].size() < s) // make sure nodes in each sCC at least larger than k = each node has at least k neighbors
            continue;

        L.push_back(make_pair(res[i], avg_vec));
        output_res.push_back(make_pair(res[i], avg_vec));
        cout << "res size: " << res[i].size() << "\tres avg: " << avg_vec << endl;

        if (avg_vec > MIC.second){
            MIC.first = res[i];
            MIC.second = avg_vec;
        }
    }
    cout << "=====> MIC size: " << MIC.first.size() << "\tinf_avg: " << MIC.second << endl;
    cout << "Init number of connect components of maximal D-core: " << output_res.size() << endl;

    //sort(L.begin(), L.end(), comparePairs);
    //// store all the nodes in the maximum d-core (d-core)
    exist_weight.resize(L.size());
    for (unsigned int i = 0; i < L.size(); i++) {
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            exist_weight[i].push_back(make_pair(L[i].first[j], weight[L[i].first[j]]));
            d_core_vec.push_back(L[i].first[j]);
        }
    }

    cout << "=====> Step3: Decomposition begin ..." << endl;
    for (unsigned int i = 0; i < L.size(); i++) {
        cout << endl << "current MIC size: " << MIC.first.size() << "\tinf_avg: " << MIC.second << endl;
        cout <<"L size is: "<< L.size() << "\tcheck L[" << i << "] now\t L[i] size is: "<< L[i].first.size()<< "\tweight: " << L[i].second << endl;
        //getchar();

        add_vec.assign(L[i].first.begin(), L[i].first.end());
        set<int> add_set(add_vec.begin(), add_vec.end());

        /*vector <int> currcore_vec;
        if (i == 0) {
            currcore_vec = d_core_vec;
        } else{
            currcore_vec = L[i-1].first;
        }
        set<int> d_core_set(currcore_vec.begin(),currcore_vec.end());
        //set<int> d_core_set(d_core_vec.begin(), d_core_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(d_core_set.begin(), d_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
         */
        //batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        //cout << "del_nodes_set.size: "<< del_nodes_set.size()<<endl;

        //////////////////////////// greedy-min-weight ////////////////////////////
        sort(exist_weight[i].begin(), exist_weight[i].end(),
             [](const pair<int, double> & left, const pair<int, double> & right){
                 return left.second < right.second;
             });
        batch_iter_del_edges_vec.clear();

        //// delete min_weight node in L[i]
        ////每次只删weight最小的点
        delete_node = exist_weight[i][0].first;
        //cout << "delete node: " << delete_node << "\t weight: " << weight[delete_node] << "\t";
        single_del_edges_vec = single_del_edges(delete_node); // delete edge L\ {L[i] \cup delete_node}
        add_set.erase(delete_node);

        unordered_set<int> rem_set(add_set.begin(), add_set.end());
        single_kl_del_edges_vec = compute_kl_core(add_set,rem_set, k, l, s);
        vector<vector<int>> resc;
        if (rem_set.size() >= s && rem_set.size() >= (k+l)){
            resc = compute_connected_components(rem_set, k, l, s);
        }
        for(const auto& component : resc) {
            if (containsVector(L, component) == false){
                /*cout << "component.size: "<< component.size() << "\tavg_weight: "<< avg_vector(component) << "\t{ ";
                for (auto it = component.begin(); it != component.end(); ++it) {
                    int node = *it;
                    cout << node <<  ' ';
                    //cout << node << '-' << weight[node] << ' ';
                }
                cout << "}"<<endl;*/
                update_weight = avg_vector(component);
                L.push_back(make_pair(component, update_weight));
                vector<pair<int, double>> nodeInfo;
                for (auto it = component.begin(); it != component.end(); ++it) {
                    int node = *it;
                    nodeInfo.push_back(make_pair(node, weight[node]));
                }
                exist_weight.push_back(nodeInfo);


                if (update_weight > MIC.second){ // keep the max weight
                    MIC.first = component;
                    MIC.second = update_weight;
                }
            }
        }

        //batch_add_edges(single_del_edges_vec); //A i=2
        //batch_add_edges(single_kl_del_edges_vec); //B A+B 对应L[1]-L[2],对应d_core_set(L[1].first.begin(), L[1].first.end());不用加了再删
        //batch_add_edges(batch_del_edges_vec); // 对应L[0]-L[1]

    }
    endTime = omp_get_wtime();
    cout <<endl << "Improved greedy avg global MIC end" << endl;

    // output result
    cout << "=====> Step4: MIC avg_inf is: "<< MIC.second << " size: " << MIC.first.size() << "\t{ ";
    for (int i = 0; i < MIC.first.size(); i ++){
        cout << MIC.first[i] << " ";
    }cout << "}"<< endl;

    cout << "k: " << k << " l: " << l << " s: " << s << endl;
    cout << "Total Time : " << (double)(endTime - startTime)  << " s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;


    // write result to file
    cout << "output..."<< endl;
    string k_name = "k_" + std::to_string(k) + "_";
    string l_name = "l_" + std::to_string(l) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    string dataset_name = get_dataset(dataset);
    //string filename = "../../result/"+dataset_name+"avg_global_mic_" +seed_mode+"_" +snum+"_" + k_name + l_name + s_name + ".txt";
    string filename = "../../result/"+dataset_name+"avg_global_mic_" +seed_mode+"_"  + p_name +snum+"_" + k_name + l_name + s_name + ".txt";
    cout << filename <<endl;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime)  << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime)  << endl;
    outfile << "Inf: "<< MIC.second << "\nSize: " << MIC.first.size() <<endl;
    for (int i = 0; i < MIC.first.size(); i ++){
        outfile << MIC.first[i] << " ";
    }
    outfile << endl;
    outfile.close();
    cout << "output end"<< endl;
}

////+++++++++++++++++++++++++++++++++++++++++++++++ local +++++++++++++++++++++++++++++++++++++++++++++++
void SMIC::improved_maxweight_avg_local_mic(const string dataset, int k, int l,int s, string seed_mode, string snum, const int p){
    cout << "====================== Improved greedy avg local MIC begin ... ======================" << endl;
    double startTime, endTime;
    double coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> L;
    vector<int> add_vec; // after delete one node, new vector
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    double update_weight; // updated weight
    int delete_node;
    double avg_vec;
    vector<int> d_core_vec; // store all the nodes in the sccs
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;

    vector<pair<int, double>> nodeset_sort;
    vector<int> node_L_loc (size_n, -1);
    pair<vector<int>, double> MIC_v;

    // init一次就可以，不需要每次执行 compute_local_max_comm_v 都init copy数组
    init_Adjlist_copy();
    init_loc_copy();

    startTime = omp_get_wtime();
    cout << "=====> Step1: Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    coreStartTime = omp_get_wtime();
    compute_max_kl_core(k,l); // compute the maximum (k,l)-core in G
    coreEndTime = omp_get_wtime();

    cout << "=====> Step2: Compute strong connected components ..." <<endl;
    res = compute_connected_components(core_nodeset, k, l, s); // compute all strong connected components (scc) on the maximum (k,l)-core

    //// construct initial L0
    MIC.second = 0.0;
    for (unsigned int i = 0; i < res.size(); i++) {
        avg_vec = avg_vector(res[i]); // compute the average weight of each sCC
        if (res[i].size() < (k+l) || res[i].size() < s) // make sure nodes in each sCC at least larger than k = each node has at least k neighbors
            continue;

        L.push_back(make_pair(res[i], avg_vec));
        for (int j = 0; j < res[i].size(); j ++){
            int node = res[i][j];
            node_L_loc[node] = L.size() - 1;
            d_core_vec.push_back(node);
        }

        cout << "res size: " << res[i].size() << "\tres avg: " << avg_vec << endl;

        if (avg_vec > MIC.second){
            MIC.first = res[i];
            MIC.second = avg_vec;
        }
    }
    cout << "=====> MIC size: " << MIC.first.size() << "\tinf_avg: " << MIC.second << endl;
    cout << "Init number of connect components of maximal D-core: " << L.size() << endl;

    cout << "=====> Step3: Sort node in compute_max_kl_core in decreasing order of weight ..." << endl;
    for (int node: core_nodeset){
        nodeset_sort.push_back(make_pair(node,weight[node]));
    }
    // sort in descreasing order of weight
    sort(nodeset_sort.begin(), nodeset_sort.end(), compare_sort);

    double epsilon = 0.1;
    cout << "=====> Step4: Decomposition begin ..." << endl;
    for (int i = 0; i < nodeset_sort.size(); i ++){
        int max_node = nodeset_sort[i].first;
        int max_node_Lid = node_L_loc[max_node];
        if (max_node_Lid == -1) // max_node_Lid == -1 <==> AdjListIn[max_node] < k || AdjListOut[max_node] < l
            continue;
        //cout << "check node: " << max_node << "\t weight: "<< weight[max_node] << "\t L id is: "<< node_L_loc[max_node]<< endl; //getchar();
        if (weight[max_node] <= MIC.second * (1 + epsilon)){
            break; //若v的weight小于当前最大avg，由于在所有包含v的社区中，v的weight是最大的，因此合并其他weight更小的节点只会使avg减小
        }
        pair<vector<int>, double> max_node_L;
        max_node_L.first = L[max_node_Lid].first;
        max_node_L.second = L[max_node_Lid].second;
        if (max_node_L.first.size() < s || max_node_L.first.size()  < (k+l))
            continue;

        ////!!!!需要删掉除max_node_L以外的其他所有子图
        set<int> Lv_set (max_node_L.first.begin(), max_node_L.first.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set<int> d_core_set (d_core_vec.begin(), d_core_vec.end()); //// d_core_vec要更新
        set_difference(d_core_set.begin(), d_core_set.end(),
                       Lv_set.begin(), Lv_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
        batch_del_edges_vec = batch_del_edges(del_nodes_vec); // delete edge L\Lv
        ////////////////// 找到包含 max_node 的所有社区中 avg influenced expectation 最大的社区 //////////////////
        cout << "=====> Step4.1: Find the maximum weight comm including node "<< max_node <<"\tweight: "<< weight[max_node] <<endl; //<<"\tindeg: "<< AdjListIn[max_node].size()<<"\toutdeg: "<< AdjListOut[max_node].size()<< endl;
        //cout << "current MIC size: " << MIC.first.size() << "\tinf_avg: " << MIC.second << endl;
        MIC_v = compute_local_max_comm_v(max_node_L, max_node, k, l, s);
        clear_copy(); //每次执行完之后清空AdjListcopy中的元素，每个内部向量都将变为空，但向量本身的大小不变

        //cout << "=====> Step4.2: MIC_v size: "<< MIC_v.first.size() << "\tinf_avg: " << MIC_v.second << endl;
        if (MIC_v.second > MIC.second){
            MIC.first = MIC_v.first;
            MIC.second = MIC_v.second;
        }
        //cout << "=====> MIC size: " << MIC.first.size() << "\tinf_avg: " << MIC.second << endl;


        ////!!! 删除当前max_node，分解原本的L[max_node_Lid]（不对，应该是分解G/v），生成新的D-core scc, 并更新其余节点的max_node_Lid
        cout << "=====> Step4.2: Delete current max_node "<< max_node << endl;
        //// 先把删掉的L\Lv 加回来
        batch_add_edges(batch_del_edges_vec);

        delete_node = max_node;
        node_L_loc[max_node] = -1;
        //cout << "delete node: " << delete_node << "\t";
        single_del_edges_vec = single_del_edges(delete_node); // delete edge L\ {L[i] \cup delete_node}


        set<int>add_set(d_core_vec.begin(), d_core_vec.end()); ///////////// 应该是在G/v上
        add_set.erase(delete_node);
        unordered_set<int> rem_set(add_set.begin(), add_set.end());
        single_kl_del_edges_vec = compute_kl_core(add_set,rem_set, k, l, s);
        //// d_core_vec要更新
        d_core_vec.assign(rem_set.begin(), rem_set.end());

        vector<vector<int>> resc;
        if (rem_set.size() >= s && rem_set.size() >= (k+l)){
            resc = compute_connected_components(rem_set, k, l, s);
        }
        set<int> ori_set(add_set.begin(), add_set.end());
        set<int> fin_set;
        set<int> del_set;

        L.clear();
        for(const auto& component : resc) {
            if (containsVector(L, component) == false){
                 //cout << "component.size: "<< component.size() << "\tavg_weight: "<< avg_vector(component) << "\t{ ";
                /*for(auto node : component) {
                    cout << node <<  ' ';
                    // cout << node << '-'<< weight[node] << ' ';
                }cout << "}"<<endl;*/
                update_weight = avg_vector(component);
                L.push_back(make_pair(component,update_weight));
                //for (int node: component) {
                for (auto it = component.begin(); it != component.end(); ++it) {
                    int node = *it;
                    node_L_loc[node] = L.size() - 1; //update location L of remain nodes
                    fin_set.insert(node);
                }
            }
        }
        set_difference(ori_set.begin(), ori_set.end(),
                       fin_set.begin(), fin_set.end(),
                       inserter(del_set, del_set.begin()));
        //for (int node : del_set){
        for (auto it = del_set.begin(); it != del_set.end(); ++it) { // 把delete_node 删除之后，不满足(k,l)条件的节点
            int node = *it;
            node_L_loc[node] = -1;
        }


    }

    cout <<endl << "Improved greedy avg local MIC end" << endl;
    endTime = omp_get_wtime();
    // output result
    cout << "=====> Step5: MIC avg_inf is: "<< MIC.second << " size: " << MIC.first.size() << "\t{ ";
    for (int i = 0; i < MIC.first.size(); i ++){
        cout << MIC.first[i] << " ";
    }
    cout << "}"<< endl;


    cout << "k: " << k << " l: " << l << " s: " << s << endl;
    //cout << "startTime: " << startTime << "\t endTime: " << endTime<< endl;
    cout << "Total Time : " << (double)(endTime - startTime)  << " s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;

    // write result to file
    cout << "output..."<< endl;
    string k_name = "k_" + std::to_string(k) + "_";
    string l_name = "l_" + std::to_string(l) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    string p_name = "p_" + std::to_string(p) + "_";
    string dataset_name = get_dataset(dataset);
    //string filename = "../../result/"+dataset_name+"avg_local_mic_" +seed_mode+"_" +snum+"_" + k_name + l_name + s_name + ".txt";
    string filename = "../../result/"+dataset_name+"avg_local_mic_" +seed_mode+"_"  + p_name +snum+"_" + k_name + l_name + s_name + ".txt";
    cout << filename <<endl;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime)  << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime)  << endl;
    outfile << "Inf: "<< MIC.second << "\nSize: " << MIC.first.size() <<endl;
    for (int i = 0; i < MIC.first.size(); i ++){
        outfile << MIC.first[i] << " ";
    }
    outfile << endl;
    outfile.close();
    cout << "output end"<< endl;

}


//// ======== min algorithm =============
void SMIC::improved_minweight_min_global_mic(const string dataset, int k, int l,int s, string seed_mode, string snum){
    cout << "====================== Improved greedy min global MIC begin ... ======================" << endl;
    double startTime, endTime;
    double coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> L; //int: connetcted component double: avg weight
    vector<pair<vector<int>, double>> res_L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> d_core_vec; // store all the nodes in the sccs
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_iter_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    vector<double> weight_vec;
    vector<vector<pair<int, double>>> exist_weight;
    double update_weight; // updated weight
    int delete_node;
    double min_vec;

    startTime = omp_get_wtime();

    cout << "=====> Step1: Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    coreStartTime = omp_get_wtime();
    compute_max_kl_core(k,l); // compute the maximum (k,l)-core in G
    coreEndTime = omp_get_wtime();

    cout << "=====> Step2: Compute strong connected components ..." <<endl;
    res = compute_connected_components(core_nodeset, k, l, s); // compute all strong connected components (scc) on the maximum (k,l)-core

    //// construct initial L0
    MIC.second = 0.0;
    for (unsigned int i = 0; i < res.size(); i++) {
        min_vec = min_vector(res[i]); // compute the min weight of each sCC
        if (res[i].size() < (k+l) || res[i].size() < s) // make sure nodes in each sCC at least larger than k = each node has at least k neighbors
            continue;

        L.push_back(make_pair(res[i], min_vec));
        output_res.push_back(make_pair(res[i], min_vec));
        cout << "res size: " << res[i].size() << "\tres min: " << min_vec << endl;

        if (min_vec >= MIC.second){
            MIC.first = res[i];
            MIC.second = min_vec;
        }
    }
    cout << "=====> MIC size: " << MIC.first.size() << "\tinf_min: " << MIC.second << endl;
    cout << "Init number of connect components of maximal D-core: " << output_res.size() << endl;

    //sort(L.begin(), L.end(), comparePairs);
    //// store all the nodes in the maximum d-core (d-core)
    exist_weight.resize(L.size());
    for (unsigned int i = 0; i < L.size(); i++) {
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            exist_weight[i].push_back(make_pair(L[i].first[j], weight[L[i].first[j]]));
            d_core_vec.push_back(L[i].first[j]);
        }
    }

    cout << "=====> Step3: Decomposition begin ..." << endl;
    for (unsigned int i = 0; i < L.size(); i++) {
        cout << endl << "current MIC size: " << MIC.first.size() << "\tinf_min: " << MIC.second << endl;
        cout <<"L size is: "<< L.size() << "\tcheck L[" << i << "] now\t L[i] size is: "<< L[i].first.size()<< "\tweight: " << L[i].second << endl;
        //getchar();

        add_vec.assign(L[i].first.begin(), L[i].first.end());
        set<int> add_set(add_vec.begin(), add_vec.end());

        /*vector <int> currcore_vec;
        if (i == 0) {
            currcore_vec = d_core_vec;
        } else{
            currcore_vec = L[i-1].first;
        }
        set<int> d_core_set(currcore_vec.begin(),currcore_vec.end());
        //set<int> d_core_set(d_core_vec.begin(), d_core_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(d_core_set.begin(), d_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
         */
        //batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        //cout << "del_nodes_set.size: "<< del_nodes_set.size()<<endl;

        //////////////////////////// greedy-min-weight ////////////////////////////
        sort(exist_weight[i].begin(), exist_weight[i].end(),
             [](const pair<int, double> & left, const pair<int, double> & right){
                 return left.second < right.second;
             });
        batch_iter_del_edges_vec.clear();

        //// delete min_weight node in L[i]
        ////每次只删weight最小的点
        delete_node = exist_weight[i][0].first;
        cout << "delete node: " << delete_node << "\t weight: " << weight[delete_node] << "\t";
        single_del_edges_vec = single_del_edges(delete_node); // delete edge L\ {L[i] \cup delete_node}
        add_set.erase(delete_node);

        unordered_set<int> rem_set(add_set.begin(), add_set.end());
        single_kl_del_edges_vec = compute_kl_core(add_set,rem_set, k, l, s);
        vector<vector<int>> resc;
        if (rem_set.size() >= s && rem_set.size() >= (k+l)){
            resc = compute_connected_components(rem_set, k, l, s);
        }
        for(const auto& component : resc) {
            if (containsVector(L, component) == false){
                /*cout << "component.size: "<< component.size() << "\tmin_weight: "<< min_vector(component) << "\t{ ";
                for (auto it = component.begin(); it != component.end(); ++it) {
                    int node = *it;
                    cout << node <<  ' ';
                    //cout << node << '-' << weight[node] << ' ';
                }
                cout << "}"<<endl;*/
                update_weight = min_vector(component);
                L.push_back(make_pair(component, update_weight));
                vector<pair<int, double>> nodeInfo;
                for (auto it = component.begin(); it != component.end(); ++it) {
                    int node = *it;
                    nodeInfo.push_back(make_pair(node, weight[node]));
                }
                exist_weight.push_back(nodeInfo);

                if (update_weight > MIC.second){ // keep the max weight
                    MIC.first = component;
                    MIC.second = update_weight;
                }
            }
        }

        //batch_add_edges(single_del_edges_vec); //A i=2
        //batch_add_edges(single_kl_del_edges_vec); //B A+B 对应L[1]-L[2],对应d_core_set(L[1].first.begin(), L[1].first.end());不用加了再删
        //batch_add_edges(batch_del_edges_vec); // 对应L[0]-L[1]

    }
    endTime = omp_get_wtime();
    cout <<endl << "Improved greedy min global MIC end" << endl;

    // output result
    cout << "=====> Step4: MIC min_inf is: "<< MIC.second << " size: " << MIC.first.size() << "\t{ ";
    for (int i = 0; i < MIC.first.size(); i ++){
        cout << MIC.first[i] << " ";
    }cout << "}"<< endl;


    cout << "k: " << k << " l: " << l << " s: " << s << endl;
    cout << "Total Time : " << (double)(endTime - startTime)  << " s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;

    // write result to file
    cout << "output..."<< endl;
    string k_name = "k_" + std::to_string(k) + "_";
    string l_name = "l_" + std::to_string(l) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    string dataset_name = get_dataset(dataset);
    string filename = "../../result/"+dataset_name+"min_global_mic_" +seed_mode+"_" +snum+"_" + k_name + l_name + s_name + ".txt";
    cout << filename <<endl;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime)  << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime)  << endl;
    outfile << "Inf: "<< MIC.second << "\nSize: " << MIC.first.size() <<endl;
    for (int i = 0; i < MIC.first.size(); i ++){
        outfile << MIC.first[i] << " ";
    }
    outfile << endl;
    outfile.close();
    cout << "output end"<< endl;

}

void SMIC::improved_one_maxweight_min_local_mic(const string dataset, int k, int l,int s, string seed_mode, string snum){
    cout << "====================== Improved greedy min local MIC begin ... ======================" << endl;
    double startTime, endTime;
    double coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> L;
    vector<int> add_vec; // after delete one node, new vector
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    double update_weight; // updated weight
    int delete_node;
    double min_vec;

    vector<pair<int, double>> nodeset_sort;
    vector<pair<int, double>> core_sort;
    vector<int> node_L_loc (size_n, -1);
    pair<vector<int>, double> MIC_v;

    // init一次就可以，不需要每次执行 compute_local_max_comm_v 都init copy数组
    init_Adjlist_copy();
    init_loc_copy();

    startTime = omp_get_wtime();
    cout << "=====> Step1: Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    coreStartTime = omp_get_wtime();
    compute_max_kl_core(k,l); // compute the maximum (k,l)-core in G
    coreEndTime = omp_get_wtime();

    cout << "=====> Step2: Compute strong connected components ..." <<endl;
    res = compute_connected_components(core_nodeset, k, l, s); // compute all strong connected components (scc) on the maximum (k,l)-core

    //// construct initial L0
    MIC.second = 0.0;
    for (unsigned int i = 0; i < res.size(); i++) {
        min_vec = min_vector(res[i]); // compute the average weight of each sCC
        if (res[i].size() < (k+l) || res[i].size() < s) // make sure nodes in each sCC at least larger than k = each node has at least k neighbors
            continue;

        L.push_back(make_pair(res[i], min_vec));

//        for (int j = 0; j < res[i].size(); j ++){
//            int node = res[i][j];
//            node_L_loc[node] = L.size() - 1;
//        }

        cout << "res size: " << res[i].size() << "\tres min: " << min_vec << endl;

        if (min_vec >= MIC.second){
            MIC.first = res[i];
            MIC.second = min_vec;
        }
    }
    cout << "=====> MIC size: " << MIC.first.size() << "\tinf_min: " << MIC.second << endl;
    cout << "Init number of connect components of maximal D-core: " << L.size() << endl;

    cout << "=====> Step3: Sort connected components in descending order of min_vec..." << endl;

    for (int node: nodeset){
        nodeset_sort.push_back(make_pair(node,weight[node]));
    }
    sort(nodeset_sort.begin(), nodeset_sort.end(), compare_sort);// sort in descreasing order of weight

    pair<vector<int>, double> max_node_L;
    max_node_L= MIC; //从当前min_inf最大的连通分量里面local的找最大的连通分量
//    for (int i = 0; i < nodeset_sort.size(); i ++){
//        int max_node = nodeset_sort[i].first;
//        int max_node_Lid = node_L_loc[max_node];
//        if (max_node_Lid == -1) // max_node_Lid == -1 <==> AdjListIn[max_node] < k || AdjListOut[max_node] < l
//            continue;
//        max_node_L.first = L[max_node_Lid].first;
//        max_node_L.second = L[max_node_Lid].second;
//
//    }

    ////////////////// 按照weight从大到小的顺序加入节点，直到找到当前 min_inf 最大的core //////////////////
    cout << "=====> Step4: Insertion begin ..." << endl;
    for (int node: max_node_L.first){
        core_sort.push_back(make_pair(node,weight[node]));
    }
    sort(core_sort.begin(), core_sort.end(), compare_sort);// sort in descreasing order of weight

    MIC_v = compute_local_max_min_comm(max_node_L, core_sort, k, l, s);
    if (MIC_v.second > MIC.second){
        MIC.first=MIC_v.first;
        MIC.second=MIC_v.second;
    }
    cout << "=====> MIC size: " << MIC.first.size() << "\tinf_min: " << MIC.second << endl;

    cout <<endl << "Improved greedy min local MIC end" << endl;

    // output result
    cout << "=====> Step5: MIC min_inf is: "<< MIC.second << " size: " << MIC.first.size() << "\t{ ";
    for (int i = 0; i < MIC.first.size(); i ++){
        cout << MIC.first[i] << " ";
    }
    cout << "}"<< endl;

    endTime = omp_get_wtime();
    cout << "k: " << k << " l: " << l << " s: " << s << endl;
    //cout << "startTime: " << startTime << "\t endTime: " << endTime<< endl;
    cout << "Total Time : " << (double)(endTime - startTime)  << " s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;

    // write result to file
    cout << "output..."<< endl;
    string k_name = "k_" + std::to_string(k) + "_";
    string l_name = "l_" + std::to_string(l) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    string dataset_name = get_dataset(dataset);
    string filename = "../../result/"+dataset_name+"min_local_mic_" +seed_mode+"_" +snum+"_" + k_name + l_name + s_name + ".txt";
    cout << filename <<endl;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime)  << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime)  << endl;
    outfile << "Inf: "<< MIC.second << "\nSize: " << MIC.first.size() <<endl;
    for (int i = 0; i < MIC.first.size(); i ++){
        outfile << MIC.first[i] << " ";
    }
    outfile << endl;
    outfile.close();
    cout << "output end"<< endl;

}

//// ======== max algorithm =============
void SMIC::improved_minweight_max_global_mic(const string dataset, int k, int l,int s, string seed_mode, string snum){
    cout << "====================== Improved greedy max global MIC begin ... ======================" << endl;
    double startTime, endTime;
    double coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> L; //int: connetcted component double: avg weight
    vector<pair<vector<int>, double>> res_L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> d_core_vec; // store all the nodes in the sccs
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_iter_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    vector<double> weight_vec;
    vector<vector<pair<int, double>>> exist_weight;
    double update_weight; // updated weight
    int delete_node;
    double max_vec;

    startTime = omp_get_wtime();

    cout << "=====> Step1: Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    coreStartTime = omp_get_wtime();
    compute_max_kl_core(k,l); // compute the maximum (k,l)-core in G
    coreEndTime = omp_get_wtime();

    cout << "=====> Step2: Compute strong connected components ..." <<endl;
    res = compute_connected_components(core_nodeset, k, l, s); // compute all strong connected components (scc) on the maximum (k,l)-core

    //// construct initial L0
    MIC.second = 0.0;
    for (unsigned int i = 0; i < res.size(); i++) {
        max_vec = max_vector(res[i]); // compute the max weight of each sCC
        if (res[i].size() < (k+l) || res[i].size() < s) // make sure nodes in each sCC at least larger than k = each node has at least k neighbors
            continue;

        L.push_back(make_pair(res[i], max_vec));
        output_res.push_back(make_pair(res[i], max_vec));
        cout << "res size: " << res[i].size() << "\tres max: " << max_vec << endl;

        if (max_vec > MIC.second){
            MIC.first = res[i];
            MIC.second = max_vec;
        }
    }
    cout << "=====> MIC size: " << MIC.first.size() << "\tinf_max: " << MIC.second << endl;
    cout << "Init number of connect components of maximal D-core: " << output_res.size() << endl;

    //sort(L.begin(), L.end(), comparePairs);
    //// store all the nodes in the maximum d-core (d-core)
    exist_weight.resize(L.size());
    for (unsigned int i = 0; i < L.size(); i++) {
        for (unsigned int j = 0; j < L[i].first.size(); j++) {
            exist_weight[i].push_back(make_pair(L[i].first[j], weight[L[i].first[j]]));
            d_core_vec.push_back(L[i].first[j]);
        }
    }

    cout << "=====> Step3: Decomposition begin ..." << endl;
    for (unsigned int i = 0; i < L.size(); i++) {
        cout << endl << "current MIC size: " << MIC.first.size() << "\tinf_max: " << MIC.second << endl;
        cout <<"L size is: "<< L.size() << "\tcheck L[" << i << "] now\t L[i] size is: "<< L[i].first.size()<< "\tweight: " << L[i].second << endl;
        //getchar();

        add_vec.assign(L[i].first.begin(), L[i].first.end());
        set<int> add_set(add_vec.begin(), add_vec.end());

        /*vector <int> currcore_vec;
        if (i == 0) {
            currcore_vec = d_core_vec;
        } else{
            currcore_vec = L[i-1].first;
        }
        set<int> d_core_set(currcore_vec.begin(),currcore_vec.end());
        //set<int> d_core_set(d_core_vec.begin(), d_core_vec.end());
        set<int> del_nodes_set;
        del_nodes_vec.clear();
        set_difference(d_core_set.begin(), d_core_set.end(),
                       add_set.begin(), add_set.end(),
                       inserter(del_nodes_set, del_nodes_set.begin()));
        del_nodes_vec.assign(del_nodes_set.begin(), del_nodes_set.end());
         */
        //batch_del_edges_vec = batch_del_edges(del_nodes_vec);
        //cout << "del_nodes_set.size: "<< del_nodes_set.size()<<endl;

        //////////////////////////// greedy-max-weight ////////////////////////////
        sort(exist_weight[i].begin(), exist_weight[i].end(),
             [](const pair<int, double> & left, const pair<int, double> & right){
                 return left.second < right.second;
             });
        batch_iter_del_edges_vec.clear();

        //// delete max_weight node in L[i]
        ////每次只删weight最小的点
        delete_node = exist_weight[i][0].first;
        cout << "delete node: " << delete_node << "\t weight: " << weight[delete_node] << "\t";
        single_del_edges_vec = single_del_edges(delete_node); // delete edge L\ {L[i] \cup delete_node}
        add_set.erase(delete_node);

        unordered_set<int> rem_set(add_set.begin(), add_set.end());
        single_kl_del_edges_vec = compute_kl_core(add_set,rem_set, k, l, s);
        vector<vector<int>> resc;
        if (rem_set.size() >= s && rem_set.size() >= (k+l)){
            resc = compute_connected_components(rem_set, k, l, s);
        }
        for(const auto& component : resc) {
            if (containsVector(L, component) == false){
                /*cout << "component.size: "<< component.size() << "\tmax_weight: "<< max_vector(component) << "\t{ ";
                for (auto it = component.begin(); it != component.end(); ++it) {
                    int node = *it;
                    cout << node <<  ' ';
                    //cout << node << '-' << weight[node] << ' ';
                }
                cout << "}"<<endl;*/
                update_weight = max_vector(component);
                L.push_back(make_pair(component, update_weight));
                vector<pair<int, double>> nodeInfo;
                for (auto it = component.begin(); it != component.end(); ++it) {
                    int node = *it;
                    nodeInfo.push_back(make_pair(node, weight[node]));
                }
                exist_weight.push_back(nodeInfo);

                if (update_weight > MIC.second){ // keep the max weight
                    MIC.first = component;
                    MIC.second = update_weight;
                }
            }
        }

        //batch_add_edges(single_del_edges_vec); //A i=2
        //batch_add_edges(single_kl_del_edges_vec); //B A+B 对应L[1]-L[2],对应d_core_set(L[1].first.begin(), L[1].first.end());不用加了再删
        //batch_add_edges(batch_del_edges_vec); // 对应L[0]-L[1]

    }
    endTime = omp_get_wtime();
    cout <<endl << "Improved greedy max global MIC end" << endl;

    // output result
    cout << "=====> Step4: MIC max_inf is: "<< MIC.second << " size: " << MIC.first.size() << "\t{ ";
    for (int i = 0; i < MIC.first.size(); i ++){
        cout << MIC.first[i] << " ";
    }cout << "}"<< endl;


    cout << "k: " << k << " l: " << l << " s: " << s << endl;
    cout << "Total Time : " << (double)(endTime - startTime)  << " s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;

    // write result to file
    cout << "output..."<< endl;
    string k_name = "k_" + std::to_string(k) + "_";
    string l_name = "l_" + std::to_string(l) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    string dataset_name = get_dataset(dataset);
    string filename = "../../result/"+dataset_name+"max_global_mic_" +seed_mode+"_" +snum+"_" + k_name + l_name + s_name + ".txt";
    cout << filename <<endl;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime)  << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime)  << endl;
    outfile << "Inf: "<< MIC.second << "\nSize: " << MIC.first.size() <<endl;
    for (int i = 0; i < MIC.first.size(); i ++){
        outfile << MIC.first[i] << " ";
    }
    outfile << endl;
    outfile.close();
    cout << "output end"<< endl;

}

void SMIC::improved_maxweight_max_local_mic(const string dataset, int k, int l,int s, string seed_mode, string snum){
    cout << "====================== Improved greedy max local MIC begin ... ======================" << endl;
    double startTime, endTime;
    double coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> L;
    vector<int> add_vec; // after delete one node, new vector
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    double update_weight; // updated weight
    int delete_node;
    double max_vec;

    vector<pair<int, double>> nodeset_sort;
    vector<int> node_L_loc (size_n, -1);
    pair<vector<int>, double> MIC_v;

    // init一次就可以，不需要每次执行 compute_local_max_comm_v 都init copy数组
    init_Adjlist_copy();
    init_loc_copy();

    startTime = omp_get_wtime();
    cout << "=====> Step1: Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    coreStartTime = omp_get_wtime();
    compute_max_kl_core(k,l); // compute the maximum (k,l)-core in G
    coreEndTime = omp_get_wtime();

    cout << "=====> Step2: Compute strong connected components ..." <<endl;
    res = compute_connected_components(core_nodeset, k, l, s); // compute all strong connected components (scc) on the maximum (k,l)-core

    //// construct initial L0
    MIC.second = 0.0;
    for (unsigned int i = 0; i < res.size(); i++) {
        max_vec = max_vector(res[i]); // compute the average weight of each sCC
        if (res[i].size() < (k+l) || res[i].size() < s) // make sure nodes in each sCC at least larger than k = each node has at least k neighbors
            continue;

        L.push_back(make_pair(res[i], max_vec));
        for (int j = 0; j < res[i].size(); j ++){
            int node = res[i][j];
            node_L_loc[node] = L.size() - 1;
        }

        cout << "res size: " << res[i].size() << "\tres max: " << max_vec << endl;

        if (max_vec > MIC.second){
            MIC.first = res[i];
            MIC.second = max_vec;
        }
    }
    cout << "=====> MIC size: " << MIC.first.size() << "\tinf_max: " << MIC.second << endl;
    cout << "Init number of connect components of maximal D-core: " << L.size() << endl;

    cout << "=====> Step3: Sort node in V in decreasing order of weight ..." << endl;
    for (int node: nodeset){
        nodeset_sort.push_back(make_pair(node,weight[node]));
    }
    // sort in descreasing order of weight
    sort(nodeset_sort.begin(), nodeset_sort.end(), compare_sort);
    /*sort(nodeset_sort.begin(), nodeset_sort.end(), [this](const pair<int, double>& a, const pair<int, double>& b) {
        if (a.second != b.second) {
            return a.second > b.second;
        } else {
            int deg_a = AdjListIn[a.first].size() + AdjListOut[a.first].size();
            int deg_b = AdjListIn[b.first].size() + AdjListOut[b.first].size();
            return deg_a > deg_b;
        }
    });*/

    double epsilon = 0.0;
    cout << "=====> Step4: Decomposition begin ..." << endl;
    for (int i = 0; i < nodeset_sort.size(); i ++){
        int max_node = nodeset_sort[i].first;
        int max_node_Lid = node_L_loc[max_node];
        if (max_node_Lid == -1) // max_node_Lid == -1 <==> AdjListIn[max_node] < k || AdjListOut[max_node] < l
            continue;
        //cout << "check node: " << max_node << "\t weight: "<< weight[max_node] << "\t L id is: "<< node_L_loc[max_node]<< endl; //getchar();
        if (weight[max_node] < MIC.second * (1 + epsilon)){
            break; //若v的weight小于当前最大max，由于在所有包含v的社区中，v的weight是最大的，因此合并其他weight更小的节点只会使max减小
        }
        pair<vector<int>, double> max_node_L;
        max_node_L.first = L[max_node_Lid].first;
        max_node_L.second = L[max_node_Lid].second;
        if (max_node_L.first.size() < s || max_node_L.first.size()  < (k+l))
            continue;


        ////////////////// 找到包含 max_node 的所有社区中 max influenced expectation 最大的社区 //////////////////
        cout << "=====> Step4.1: Find the maximum weight comm including node "<< max_node <<"\tweight: "<< weight[max_node] <<endl; //<<"\tindeg: "<< AdjListIn[max_node].size()<<"\toutdeg: "<< AdjListOut[max_node].size()<< endl;
        //cout << "current MIC size: " << MIC.first.size() << "\tinf_max: " << MIC.second << endl;
        MIC_v = compute_local_max_comm_v(max_node_L, max_node, k, l, s);
        clear_copy(); //每次执行完之后清空AdjListcopy中的元素，每个内部向量都将变为空，但向量本身的大小不变

        cout << "=====> Step4.2: MIC_v size: "<< MIC_v.first.size() << "\tinf_max: " << MIC_v.second << endl;
        if (MIC_v.second > MIC.second){
            MIC.first = MIC_v.first;
            MIC.second = MIC_v.second;
        }
        cout << "=====> MIC size: " << MIC.first.size() << "\tinf_max: " << MIC.second << endl;


        // 删除当前max_node，分解原本的L[max_node_Lid]，生成新的D-core scc, 并更新其余节点的max_node_Lid
        cout << "=====> Step4.2: Delete current max_node "<< max_node << endl;
        delete_node = max_node;
        node_L_loc[max_node] = -1;
        //cout << "delete node: " << delete_node << "\t";
        single_del_edges_vec = single_del_edges(delete_node); // delete edge L\ {L[i] \cup delete_node}

        add_vec.assign(max_node_L.first.begin(), max_node_L.first.end());
        set<int>add_set(add_vec.begin(), add_vec.end());
        add_set.erase(delete_node);
        unordered_set<int> rem_set(add_set.begin(), add_set.end());
        single_kl_del_edges_vec = compute_kl_core(add_set,rem_set, k, l, s);
        vector<vector<int>> resc;
        if (rem_set.size() >= s && rem_set.size() >= (k+l)){
            resc = compute_connected_components(rem_set, k, l, s);
        }
        set<int> ori_set(add_vec.begin(), add_vec.end());
        set<int> fin_set;
        set<int> del_set;

        for(const auto& component : resc) {
            if (containsVector(L, component) == false){
                // cout << "component.size: "<< component.size() << "\tmax_weight: "<< max_vector(component) << "\t{ ";
                /*for(auto node : component) {
                    cout << node <<  ' ';
                    // cout << node << '-'<< weight[node] << ' ';
                }cout << "}"<<endl;*/
                update_weight = max_vector(component);
                L.push_back(make_pair(component,update_weight));
                //for (int node: component) {
                for (auto it = component.begin(); it != component.end(); ++it) {
                    int node = *it;
                    node_L_loc[node] = L.size() - 1; //update location L of remain nodes
                    fin_set.insert(node);
                }
            }
        }
        set_difference(ori_set.begin(), ori_set.end(),
                       fin_set.begin(), fin_set.end(),
                       inserter(del_set, del_set.begin()));
        //for (int node : del_set){
        for (auto it = del_set.begin(); it != del_set.end(); ++it) { // 把delete_node 删除之后，不满足(k,l)条件的节点
            int node = *it;
            node_L_loc[node] = -1;
        }


    }

    cout <<endl << "Improved greedy max local MIC end" << endl;
    endTime = omp_get_wtime();
    // output result
    cout << "=====> Step5: MIC max_inf is: "<< MIC.second << " size: " << MIC.first.size() << "\t{ ";
    for (int i = 0; i < MIC.first.size(); i ++){
        cout << MIC.first[i] << " ";
    }
    cout << "}"<< endl;


    cout << "k: " << k << " l: " << l << " s: " << s << endl;
    //cout << "startTime: " << startTime << "\t endTime: " << endTime<< endl;
    cout << "Total Time : " << (double)(endTime - startTime)  << " s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;

    // write result to file
    cout << "output..."<< endl;
    string k_name = "k_" + std::to_string(k) + "_";
    string l_name = "l_" + std::to_string(l) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    string dataset_name = get_dataset(dataset);
    string filename = "../../result/"+dataset_name+"max_local_mic_" +seed_mode+"_" +snum+"_" + k_name + l_name + s_name + ".txt";
    cout << filename <<endl;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime)  << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime)  << endl;
    outfile << "Inf: "<< MIC.second << "\nSize: " << MIC.first.size() <<endl;
    for (int i = 0; i < MIC.first.size(); i ++){
        outfile << MIC.first[i] << " ";
    }
    outfile << endl;
    outfile.close();
    cout << "output end"<< endl;

}

//// ======== sum algorithm =============
void SMIC::improved_sum_mic(const string dataset, int k, int l,int s, string seed_mode, string snum){
    cout << "====================== Improved greedy sum global MIC begin ... ======================" << endl;
    double startTime, endTime;
    double coreStartTime, coreEndTime;
    vector<vector<int>> res;
    vector<pair<vector<int>, double>> L; //int: connetcted component double: sum weight
    vector<pair<vector<int>, double>> res_L;
    vector<pair<vector<int>, double>> output_res;
    vector<int> add_vec; // after delete one node, new vector
    vector<int> d_core_vec; // store all the nodes in the sccs
    vector<int> del_nodes_vec;
    vector<pair<int, int>> batch_del_edges_vec;
    vector<pair<int, int>> batch_iter_del_edges_vec;
    vector<pair<int, int>> single_del_edges_vec;
    vector<pair<int, int>> single_kl_del_edges_vec;
    vector<double> weight_vec;
    vector<vector<pair<int, double>>> exist_weight;
    double update_weight; // updated weight
    int delete_node;
    double sum_vec;

    startTime = omp_get_wtime();

    cout << "=====> Step1: Compute max ("<<k<<","<<l<<")-core ..." <<endl;
    coreStartTime = omp_get_wtime();
    compute_max_kl_core(k,l); // compute the maximum (k,l)-core in G
    coreEndTime = omp_get_wtime();

    cout << "=====> Step2: Compute strong connected components ..." <<endl;
    res = compute_connected_components(core_nodeset, k, l, s); // compute all strong connected components (scc) on the maximum (k,l)-core

    //// construct initial L0
    MIC.second = 0.0;
    for (unsigned int i = 0; i < res.size(); i++) {
        sum_vec = sum_vector(res[i]); // compute the max weight of each sCC
        if (res[i].size() < (k+l) || res[i].size() < s) // make sure nodes in each sCC at least larger than k = each node has at least k neighbors
            continue;

        L.push_back(make_pair(res[i], sum_vec));
        output_res.push_back(make_pair(res[i], sum_vec));
        cout << "res size: " << res[i].size() << "\tres sum: " << sum_vec << endl;

        if (sum_vec > MIC.second){
            MIC.first = res[i];
            MIC.second = sum_vec;
        }
    }
    cout << "=====> MIC size: " << MIC.first.size() << "\tinf_sum: " << MIC.second << endl;
    cout << "Init number of connect components of maximal D-core: " << output_res.size() << endl;

    endTime = omp_get_wtime();
    cout <<endl << "Improved greedy sum global MIC end" << endl;

    // output result
    cout << "=====> Step4: MIC sum_inf is: "<< MIC.second << " size: " << MIC.first.size() << "\t{ ";
    for (int i = 0; i < MIC.first.size(); i ++){
        cout << MIC.first[i] << " ";
    }cout << "}"<< endl;


    cout << "k: " << k << " l: " << l << " s: " << s << endl;
    cout << "Total Time : " << (double)(endTime - startTime)  << " s" << endl;
    cout << "Core Time : " << (double)(coreEndTime - coreStartTime) << " s" << endl;

    // write result to file
    cout << "output..."<< endl;
    string k_name = "k_" + std::to_string(k) + "_";
    string l_name = "l_" + std::to_string(l) + "_";
    string s_name = "s_" + std::to_string(s) + "_";
    string dataset_name = get_dataset(dataset);
    string filename = "../../result/"+dataset_name+"sum_mic_" +seed_mode+"_" +snum+"_" + k_name + l_name + s_name + ".txt";
    cout << filename <<endl;
    ofstream outfile(filename);
    outfile << "Total Time (s): " << endl;
    outfile << (double)(endTime - startTime)  << endl;
    outfile << "Core Time (s): " << endl;
    outfile << (double)(coreEndTime - coreStartTime)  << endl;
    outfile << "Inf: "<< MIC.second << "\nSize: " << MIC.first.size() <<endl;
    for (int i = 0; i < MIC.first.size(); i ++){
        outfile << MIC.first[i] << " ";
    }
    outfile << endl;
    outfile.close();
    cout << "output end"<< endl;

}









