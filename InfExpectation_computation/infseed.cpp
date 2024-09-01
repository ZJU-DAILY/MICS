//
// Created by DBL-XQ on 2024/3/5.
//
#include "option.h"
#include "graph_inf.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <random>
#include <omp.h>

using namespace std;

typedef pair<float,int> fi;
typedef vector<fi> vfi;
const unsigned int UI_MAX = 4294967295U;

int generateInfluenceSampleLT(Graph &g,vector<bool> &visit, vector<int> &visit_index, vector<double> &threshold, int root)
{
    int i, cur;
    float flip;

    int curPos = 0;
    int num_marked = 1;
    visit[root] = true;
    visit_index[0] = root;
    threshold[root] = -1.0;

    while(curPos < num_marked) {
        cur = visit_index[curPos];
        const vector<int> &neigh = g.OutNode[cur];
        for (i = 0; i < g.outDeg[cur]; i++) {
            if(threshold[neigh[i]] < 0)
                continue;
            if(threshold[neigh[i]] > 1.0){
                flip = sfmt_genrand_uint32(&g.sfmtSeed) / (float)UI_MAX;
                threshold[neigh[i]] = flip;
            }
            else
                flip = threshold[neigh[i]];
            const vector<int> &inNeigh = g.InNode[neigh[i]];
            unsigned int inD = g.inDeg[neigh[i]];
            for(int k = 0 ; k < inD ; k++){
                int o = inNeigh[k];
                // for activated node
                if (threshold[o] < 0.0){
                    flip -= g.probT[neigh[i]][k];
                }
            }
            if (flip <= 0.0) {
                threshold[neigh[i]] = -1.0;
                if (!visit[neigh[i]]) {
                    visit_index[num_marked] = neigh[i];
                    num_marked++;
                    visit[neigh[i]] = true;
                }
            }
        }
        curPos++;
    }

    for(i = 0; i < num_marked; i++) {
        visit[visit_index[i]] = false;
    }
    for(i = 0; i < (int)threshold.size(); i++) {
        threshold[i] = 2.0;
    }
    return num_marked;
}

int generateInfluenceSampleIC(Graph &g,vector<bool> &visit, vector<int> &visit_index, int root)
{
    if(root >= g.prob.size())
        cout << "root = " << root << " enter\n";
    int i, cur;
    float flip;

    int curPos = 0;
    int num_marked = 1;
    visit[root] = true;
    visit_index[0] = root;

    while(curPos < num_marked) {
        cur = visit_index[curPos];
        if(cur >= g.prob.size())
        {
            cout <<"WRONG, cur = " << cur << endl;
        }
        for (i = 0; i < g.prob[cur].size(); i++) {//i < g.outDeg[cur]
            flip = sfmt_genrand_uint32(&g.sfmtSeed) / (float)UI_MAX;

            if(flip <= g.prob[cur][i].second &&  !visit[g.prob[cur][i].first]){ //flip <= g.prob[cur][neigh[i]]
                visit_index[num_marked] = g.prob[cur][i].first;
                num_marked++;
                visit[g.prob[cur][i].first] = true;
            }
        }
        curPos++;
    }

    for(i = 0; i < num_marked; i++) {
        visit[visit_index[i]] = false;
    }
    return num_marked;
}

bool compareDescending(const std::pair<int, double>& a, const std::pair<int, double>& b) {
    return a.second > b.second;
}

// estimate influence of singleton seeds using 10K monte carlo simulations
void estimateInfluence(Graph &g, set<int> nodeset, int n, int r, vector<pair<int, double>> &inf) {
    vector<int> nodes;
    copy(nodeset.begin(), nodeset.end(), back_inserter(nodes));
//    for(int i = 0; i < nodes.size(); i++){
//        if(nodes[i] > nodes.size()){
//            cout << nodes[i] << endl;
//        }
//    }
    cout << "nodes.size = " << nodes.size() << endl;
    cout << "n = " << n << endl;
    //cout << nodes[1632802] << " " << nodes[0] << endl;
    cout << g.prob.size() << endl;
//    omp_set_num_threads(2);
//	#pragma omp parallel
//	{
    long long int running_total;
    float influence;
    vector<bool> visit(n,false);
    vector<int> visit_index(n,0);
    vector<double> threshold(n,2.0);
//		#pragma omp for
    for (int i = 0; i < n ; i++) {
        running_total = 0;
        for (int j = 0; j < r; j++) {
            if(g.influModel == Graph::LT)
                running_total += generateInfluenceSampleLT(g, visit, visit_index, threshold, i);
            else if(g.influModel == Graph::IC){
//                    cout << "Yes\n";
                running_total += generateInfluenceSampleIC(g, visit, visit_index, i);
            }
        }
        influence = (double)running_total / r;
        inf.push_back(make_pair(i,influence));
    }
//	}
    cout << "sort inf"<<endl;
    sort(inf.begin(), inf.end(),compareDescending);
}

int main(int argc, char ** argv)
{

    if (argc == 1){
        cout << "Usage: ./InfExpectation_computation --r 10000 --dataset ../dataset/Wiki-Vote/ --model IC " << endl;
        exit(0);
    }

    string model;
    string dataset;
    int r = 10000;

    for (int i = 1; i<argc; i++){
        if (string(argv[i]) == "--r")
            r = atoi(argv[i + 1]);
        if (string(argv[i]) == "--dataset")
            dataset = argv[i + 1];
        if (string(argv[i]) == "--model")
            model = argv[i + 1];
    }

    ASSERT((model == "IC") || (model == "LT"));

    string dir,filename;

    dir = dataset;
    filename = dataset + "graph.txt";
    Graph g(dir, filename);
    if(model == "IC")
        g.influModel = Graph::IC;
    else
        g.influModel = Graph::LT;
    int n = g.size_n;


    cout << "\n*******************" << endl;
    cout << "\tSTART" << endl;
    cout << "*******************\n" << endl;

    set<int> nodeset;
    copy(g.nodes.begin(), g.nodes.end(), inserter(nodeset, nodeset.end()));

    //vector<double> inf (n, 0.0);
    vector<pair<int, double>> inf;
    double start = omp_get_wtime();

    estimateInfluence(g, nodeset, n, r, inf);
    cout << "Time to estimate singleton influences: " << (omp_get_wtime()-start) << "s" << endl;
    string outFile = dataset +"node_inf_"+ model +".txt";
    ofstream of(outFile);

    for (int i = 0; i <inf.size(); i ++){
        of << inf[i].first << "\t" <<  inf[i].second << endl;
        //cout << inf[i].first << "\t" <<  inf[i].second << endl;

    }
    of.close();

    cout << "\n*******************" << endl;
    cout << "\tALL DONE" << endl;
    cout << "*******************" << endl;

    //
    // ALL DONE
    //
}