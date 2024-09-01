//
// Created by DBL-XQ on 2024/3/18.
//

#ifndef INC_20240306_INFEXPECTATION_COMPUTATION_GRAPH_INF_H
#define INC_20240306_INFEXPECTATION_COMPUTATION_GRAPH_INF_H

#endif //INC_20240306_INFEXPECTATION_COMPUTATION_GRAPH_INF_H

#include <algorithm>
#include <iostream>
#include <cstring>
#include <cstdint>
#include <sstream>
#include <stdio.h>
// for mmap:
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

//for customer
#include "sfmt/SFMT.h"
#include <vector>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>
#include <set>
#include <queue>
#include <unistd.h>  //close open
#include <fstream>

#define ASSERT(v) {if (!(v)) {cerr<<"ASSERT FAIL @ "<<__FILE__<<":"<<__LINE__<<endl; exit(1);}}

using namespace std;

//typedef double (*pf)(int, int);

void handle_error(const char* msg);

class Graph
{
public:
    int n, m;
    int size_n = 0;
    set<int> nodes;

    vector<int> inDeg;
    vector<int> outDeg;
    vector<vector<int>> InNode;
    vector<vector<int>> OutNode;

    vector<vector<double>> probT;
    vector<vector<pair<int,double>>> prob;
    vector<vector<int>> prob_neigh;

    sfmt_t sfmtSeed;

    enum InfluModel {IC, IC_van, LT};
    InfluModel influModel;
    void setInfuModel(InfluModel p)
    {
        influModel = p;
    }


    string folder;
    string graph_file;
    string per;
    void readNM()
    {
//        cout << "Read attribute ..." << endl;
//        std::ifstream infile((folder + "attribute.txt"));
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
    }

    void add_edge(int a, int b)
    {
        InNode[b].push_back(a);
        OutNode[a].push_back(b);
        inDeg[b]++;
        outDeg[a]++;
    }

    vector<bool> hasnode;

    void readGraph()
    {
        cout << "Read graph_inf ..." << endl;
        std::ifstream infile(graph_file);
        if (!infile.is_open())
        {
            std::cout << "The file \"" + graph_file + "\" can NOT be opened\n";
            return;
        }
        for (unsigned int i = 0; i < m; i++){
            unsigned int a, b;
            infile >> a >> b;
            //cout << "a: " << a << "\tb: " <<b<< endl;
            if (a == b)
                continue;
//            hasnode[a]=true;
//            hasnode[b]=true;
            add_edge(a, b);
            if(a>=size_n) size_n=a;
            if(b>=size_n) size_n=b;
            nodes.insert(a); nodes.insert(b);
        }
        size_n ++;
        infile.close();

        //// probability on edge : Weighted cascade （1/d_v^in）
        set<int>::iterator Iter;
        cout << "max_n: " << size_n<< endl;
        cout << "nodes.size: "<<nodes.size()<<endl;
        //probT = vector<vector<double>>(size_n);
        prob = vector<vector<pair<int,double>>>(size_n);
        //prob_neigh = vector <vector <int>> (size_n);
        //probT = vector<vector<double>>(size_n+1, vector<double>(size_n+1,0.0));
        //prob = vector<vector<double>>(size_n+1, vector<double>(size_n+1,0.0));
        for(Iter = nodes.begin(); Iter != nodes.end(); Iter ++)
        {
            int Idx = *Iter;
            double weight;
            if (inDeg[Idx] == 0) {weight = 0.; continue;}
            else weight = 1.0 / inDeg[Idx];
            //cout << "node id: " << Idx  << "\tinDeg: " << inDeg[Idx] << "\tprobability weight: " << weight << endl;
            for (int i = 0; i < InNode[Idx].size(); i ++)
            {
                int Idy = InNode[Idx][i]; //Idy: in_neighbor Idx: out_neighbor
                //probT[Idx].push_back(weight);
                prob[Idy].push_back(make_pair(Idx,weight));
                //prob[Idy].push_back(weight);
                //prob_neigh[Idy].push_back(Idx);
            }
        }
    }


    Graph(string folder, string graph_file, string per): folder(folder), graph_file(graph_file), per(per)
    {
        srand(time(NULL));
        sfmt_init_gen_rand(&sfmtSeed, rand());

        readNM();

        for (int i = 0; i < n*2; i++)
        {

            OutNode.push_back(vector<int>());
            InNode.push_back(vector<int>());
//            hasnode.push_back(false);
            inDeg.push_back(0);
            outDeg.push_back(0);
        }
        readGraph();

    }

};

double sqr(double t)
{
    return t * t;
}

void handle_error(const char* msg) {
    perror(msg);
    exit(255);
}

