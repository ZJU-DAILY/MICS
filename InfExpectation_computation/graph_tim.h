//
// Created by DBL-XQ on 2024/3/19.
//

#ifndef INC_20240306_INFEXPECTATION_COMPUTATION_GRAPH_TIM_H
#define INC_20240306_INFEXPECTATION_COMPUTATION_GRAPH_TIM_H

#endif //INC_20240306_INFEXPECTATION_COMPUTATION_GRAPH_TIM_H
#define HEAD_INFO
#include "sfmt/SFMT.h"
using namespace std;
typedef double (*pf)(int,int);
void handle_error(const char* msg);

class Graph
{
public:
    int n, m, k, size_n;
    vector<int> inDeg;
    vector<vector<int>> gT;

    vector<vector<double>> probT;

    vector<vector<int>> InNode;
    set<int> nodes;
    vector<bool> visit;
    vector<int> visit_mark;
    enum InfluModel {IC, LT};
    InfluModel influModel;
    void setInfuModel(InfluModel p){
        influModel=p;
    }



    string folder;
    string graph_file;
    string per;

    void readNM()
    {
        string attr="attribute.txt";
        cout << "Read attribute ..." << endl;
        if (per == "100"){
            attr="attribute.txt";
        }
        else{
            attr="attribute_p"+per+".txt";
        }
        std::ifstream infile((folder + attr));

        if (!infile.is_open())
        {
            std::cout << "The file \"" + graph_file + "\" can NOT be opened\n";
            return;
        }
        infile >> n >> m;
        cout << "n: " << n << "\tm: " <<m<< endl;
        infile.close();
    }

    void add_edge(int a, int b, double p){
        probT[b].push_back(p);
        gT[b].push_back(a);
        inDeg[b]++;
    }

    vector<bool> hasnode;
    void readGraph()
    {
        cout << "Read graph ..." << endl;
        size_n = 0;
        std::ifstream infile(graph_file);
        if (!infile.is_open())
        {
            std::cout << "The file \"" + graph_file + "\" can NOT be opened\n";
            return;
        }
        for (unsigned int i = 0; i < m; i++){
            unsigned int a, b;
            infile >> a >> b;
            if (a == b)
                continue;
            hasnode[a]=true;
            hasnode[b]=true;
            //add_edge(a, b);
            gT[b].push_back(a);
            inDeg[b]++;
            if(a>=size_n) size_n=a;
            if(b>=size_n) size_n=b;
            nodes.insert(a); nodes.insert(b);
            //cout << "a: " << a << "\tb: " <<b<< endl;
        }
        infile.close();

        visit_mark=vector<int>(size_n);
        visit=vector<bool>(size_n);

//// probability on edge : Weighted cascade （1/d_v^in）
        set<int>::iterator Iter;
        cout << "max_n: " << size_n<< endl;
        cout << "nodes.size: "<<nodes.size()<<endl;
        for(Iter = nodes.begin(); Iter != nodes.end(); Iter ++)
        {
            int Idx = *Iter;
            double weight;
            if (inDeg[Idx] == 0) {weight = 0.; continue;}
            else weight = 1.0 / inDeg[Idx];
            // cout << "node id: " << Idx  << "\tinDeg: " << inDeg[Idx] << "\tprobability weight: " << weight << endl;
            for (int i = 0; i < gT[Idx].size(); i ++)
            {
                probT[Idx].push_back(weight); //probT[b].push_back(p);

            }
        }
    }




    Graph(string folder, string graph_file, string per):folder(folder), graph_file(graph_file), per(per){

        readNM();

        //init vector
//        FOR(i, n){
//            gT.push_back(vector<int>());
//            hasnode.push_back(false);
//            probT.push_back(vector<double>());
//            //hyperGT.push_back(vector<int>());
//            inDeg.push_back(0);
//        }
        for (int i = 0; i < n*2; i++)
        {
            InNode.push_back(vector<int>());

            gT.push_back(vector<int>());
            hasnode.push_back(false);
            probT.push_back(vector<double>());
            //hyperGT.push_back(vector<int>());
            inDeg.push_back(0);
        }

        readGraph();

        n = size_n+1; //防止越界
    }

};
double sqr(double t)
{
    return t*t;
}
void handle_error(const char* msg) {
    perror(msg);
    exit(255);
}
#include "infgraph.h"
#include "timgraph.h"