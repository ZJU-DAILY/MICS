//
// Created by DBL-XQ on 2024/2/12.
//

#ifndef INC_20240207_INFEXPECTATION_COMPUTATION_INFCOMPUTE_H
#define INC_20240207_INFEXPECTATION_COMPUTATION_INFCOMPUTE_H

#endif //INC_20240207_INFEXPECTATION_COMPUTATION_INFCOMPUTE_H

#include <cmath>
#include "graph_inf.h"
#include <omp.h>


class InfC : public Graph
{
public:
    InfC(string folder, string graph_file, string per) :Graph(folder, graph_file, per){
    }

    set<int> seedSet;

    vector<double> Influence(int theta){// if not node then edge
        vector<double> influence;
        if (influModel == IC_van){
            influence = InfluenceIC_van(theta, seedSet);
        }
        else if (influModel == IC){
            influence = InfluenceIC_sub(theta, seedSet);
        }
        else if (influModel == LT){
            influence = InfluenceLT(theta, seedSet);

        }
        else  ASSERT(false);

        return influence;
    }


    vector<double> InfluenceIC_van(int theta, set<int> startSet){
        ASSERT(theta > 0);
        vector<double> inf (size_n, 0.0);

        for (int x = 0; x<theta; x++){
            vector< int> visit_local(size_n,false);
            queue<int> q;

            for (auto it = startSet.begin(); it != startSet.end(); ++it) {
                int seed = *it;
                q.push(seed);
                visit_local[seed] = true;
                inf[seed] ++; // compute influenced expectation Ie(v,S)
            }

            /*for (int seed : startSet){
                q.push(seed);
                visit_local[seed] = true;
                inf[seed] ++; // compute influenced expectation Ie(v,S)
            }*/

            while (!q.empty()) {
                int expand = q.front();
                q.pop();
                int node = expand;
                for (int j = 0; j < (int)prob[node].size(); j ++){ //j<(int)OutNode[node].size();
                    int v = prob[node][j].first; //v = OutNode[node][j];
                    double randDouble = sfmt_genrand_real1(&sfmtSeed);
                    //cout << " randDouble: "<< randDouble<< " probT: "<< probT[node][v];

                    if (randDouble > prob[node][j].second)//randDouble > probT[node][j]
                        continue;

                    if (visit_local[v])
                        continue;

                    visit_local[v] = true;
                    inf[v] ++; // compute influenced expectation Ie(v,S)
                    q.push(v);
                }
            }
        }
        for (int node : nodes){
            inf[node] = inf[node]/theta *1.0;
        }

        return inf;
    }


    inline double Logarithm(const double x)
    {
        // log2f is 10% faster
        return log2f(x);
    }

    vector<double> InfluenceIC_sub(int theta, set<int> startSet){
        ASSERT(theta > 0);
        vector<double> inf (size_n, 0.0);
        //double time1 = 0.0, time2 = 0.0, time3 = 0.0, time4 = 0.0, time5 = 0.0;

        // 节点的out-neigh按照weight降序排列
        for (int i : nodes){
            sort(prob[i].begin(), prob[i].end(), [](const pair<int, double>& a, const pair<int, double>& b) {
                return a.second > b.second; // 按照第二个元素（double 值）降序排列
            });
        }

        double p_threshold = 0.1;

        queue<int> q;
        for (int x = 0; x<theta; x++){

            vector<bool> visit_local(size_n,false);

            for (auto it = startSet.begin(); it != startSet.end(); ++it) {
                int seed = *it;
                q.push(seed);
                visit_local[seed] = true;
                inf[seed] ++; // compute influenced expectation Ie(v,S)
            }

            while (!q.empty()) {
                int node = q.front();
                q.pop();
                //int out_degree = (int)prob[node].size();
                //if (out_degree <= 0)  continue;
                int startMin = 0;
                while (startMin < prob[node].size()) {
                    const auto &currentedge = prob[node][startMin];
                    int nbrId = currentedge.first;
                    double node_prob = currentedge.second;
                    startMin++;

                    if (visit_local[nbrId]) continue;

                    if (node_prob < p_threshold)
                    {
                        break;
                    }

                    double randDouble = sfmt_genrand_real1(&sfmtSeed);
                    //startMin++;

                    if (randDouble > node_prob) continue;
                    //if (visit_local[nbrId]) continue;

                    visit_local[nbrId] = true;
                    inf[nbrId]++; // compute influenced expectation Ie(v,S)
                    q.push(nbrId);

                }

                while (startMin < prob[node].size()) // skip  part of out-neighbors
                {
                    double bucket_probability = prob[node][startMin].second;
                    const double log_prob = Logarithm(1 - bucket_probability); //log2f(x)
                    double _prob = sfmt_genrand_real1(&sfmtSeed);
                    startMin += floor(Logarithm(_prob) / log_prob);

                    if (startMin >= prob[node].size()) {
                        break;
                    }

                    const auto &currentedge = prob[node][startMin];
                    int nbrId = currentedge.first;
                    double accept_probability = currentedge.second;
                    double randDouble = sfmt_genrand_real1(&sfmtSeed);
                    startMin++;

                    if (randDouble > accept_probability / bucket_probability) continue;
                    if (visit_local[nbrId]) continue;

                    visit_local[nbrId] = true;
                    inf[nbrId]++; // compute influenced expectation Ie(v,S)
                    q.push(nbrId);
                }

            }
        }

        for (int node : nodes){
            inf[node] = inf[node]/theta *1.0;
        }



        return inf;
    }


    vector<double> InfluenceLT(int theta, set<int> startSet){
        ASSERT(theta > 0);
        vector<double> inf (size_n);
        //vector<double> Pr (theta, 1.0);
        vector<double> active(size_n);


        for (int x = 0; x<theta; x++){

            //for (int i = 0; i<size_n; i++){
            for (auto it = nodes.begin(); it != nodes.end(); ++it) {
                int i = *it;
                double randDouble = sfmt_genrand_real1(&sfmtSeed);
                active[i] = randDouble; //the threshold is initialized every time for each Forward Influence sampling.
            }

            queue<int> q;
            vector<bool> visit_local(size_n, false);
            for (auto it = startSet.begin(); it != startSet.end(); ++it) {
                int seed = *it;
                q.push(seed);
                visit_local[seed] = true;
                inf[seed] ++; // compute influenced expectation Ie(v,S)
            }

            while (!q.empty()) {
                int expand = q.front();
                q.pop();
                for (auto it = prob[expand].begin(); it != prob[expand].end(); ++it) {
                    pair <int, double> out_neibor = *it;
                    int u = out_neibor.first; //OutNode[u][j]
                    double p = out_neibor.second;
                    if (active[u]>0 && active[u] - p<0){ // v is activated
                        if (visit_local[u] == true)
                            continue;

                        visit_local[u] = true;
                        inf[u] ++; // compute influenced expectation Ie(v,S)
                        q.push(u);

                    }
                    active[u] -= p;
                }
            }

        }
        for (int node : nodes){
            inf[node] = inf[node]/theta;
        }
        return inf;
    }


/*
    vector<double> InfluenceLT(int theta, set<int> startSet){
        ASSERT(theta > 0);
        vector<double> inf (size_n);
        //vector<double> Pr (theta, 1.0);
        vector<double> active(size_n);


        for (int x = 0; x<theta; x++){

            //for (int i = 0; i<size_n; i++){
            for (auto it = nodes.begin(); it != nodes.end(); ++it) {
                int i = *it;
                double randDouble = sfmt_genrand_real1(&sfmtSeed);
                active[i] = randDouble; //the threshold is initialized every time for each Forward Influence sampling.
            }

            queue<int> q;
            vector<bool> visit_local(size_n, false);
            for (auto it = startSet.begin(); it != startSet.end(); ++it) {
                int seed = *it;
                q.push(seed);
                visit_local[seed] = true;
                inf[seed] ++; // compute influenced expectation Ie(v,S)
            }

            while (!q.empty()) {
                int expand = q.front();
                q.pop();
                set<int> t = ExpandLT(expand, active);

                for (int u : t){
                    if (visit_local[u] == true)
                        continue;

                    visit_local[u] = true;
                    inf[u] ++; // compute influenced expectation Ie(v,S)
                    q.push(u);
                }

            }

        }
        for (int node : nodes){
            inf[node] = inf[node]/theta;
        }
        return inf;
    }

    set<int> ExpandLT(int u, vector<double> & active){
        set<int> rtn;
        //for (int j = 0; j < prob[u].size(); j++){ //j < OutNode[u].size()
        for (auto it = prob[u].begin(); it != prob[u].end(); ++it) {
            pair <int, double> out_neibor = *it;
            int v = out_neibor.first; //OutNode[u][j]
            double p = out_neibor.second;
            if (active[v]>0 && active[v] - p<0){
                rtn.insert(v);

            }
            active[v] -= p;
        }
        return rtn;
    }

*/




};

