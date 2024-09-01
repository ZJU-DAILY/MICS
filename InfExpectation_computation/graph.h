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
    vector<vector<double>> prob;

	sfmt_t sfmtSeed;

    enum InfluModel {IC, LT};
    InfluModel influModel;
    void setInfuModel(InfluModel p)
    {
        influModel = p;
    }

    string folder;
    string graph_file;
    void readNM()
    {
        cout << "Read attribute ..." << endl;
        std::ifstream infile((folder + "attribute.txt"));
        if (!infile.is_open())
        {
            std::cout << "The file \"" + graph_file + "\" can NOT be opened\n";
            return;
        }
        infile >> n >> m;
        cout << "n: " << n << "\tm: " <<m<< endl;
        infile.close();
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
        cout << "Read graph ..." << endl;
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
            add_edge(a, b);
            if(a>=size_n) size_n=a;
            if(b>=size_n) size_n=b;
            nodes.insert(a); nodes.insert(b);
            //cout << "a: " << a << "\tb: " <<b<< endl;
        }

        infile.close();

//        //// probability on edge : Weighted cascade （1/d_v^in）
//        set<int>::iterator Iter;
//        cout << "max_n: " << size_n<< endl;
//        cout << "nodes.size: "<<nodes.size()<<endl;
//        probT = vector<vector<double>>(size_n+1, vector<double>(size_n+1,0.0));
//        //prob = vector<vector<double>>(size_n+1, vector<double>(size_n+1,0.0));
//        for(Iter = nodes.begin(); Iter != nodes.end(); Iter ++)
//        {
//            int Idx = *Iter; cout << "node: "<< Idx <<endl;
//            double weight;
//            if (inDeg[Idx] == 0) {weight = 0.; continue;}
//            else weight = 1.0 / inDeg[Idx];
//            cout << "node id: " << Idx  << "\tinDeg: " << inDeg[Idx] << "\tprobability weight: " << weight << endl;
//            for (int i = 0; i < InNode[Idx].size(); i ++)
//            {
//                int Idy = InNode[Idx][i]; //Idy: in_neighbor Idx: out_neighbor
//                probT[Idy][Idx] = weight;
//                //prob[Idx][Idy] = weight;
//            }
//        }
/*
        ofstream out("../../dataset/dblp/graph_prob.txt");
        set <int> :: iterator iter;
        int max_deg = 0;
        for(iter = nodes.begin(); iter != nodes.end(); iter ++){
            auto node = *iter;
            if (gT[node].size()==0) continue;
            else{
                for (int j = 0; j < gT[node].size(); j++){
                    out << gT[node][j] << " " << node <<" "<< probT[node][j] << endl;
                    cout << gT[node][j] << " " << node <<" "<< probT[node][j] << endl;
                }
            }
        }
        out.close();*/
    }


    Graph(string folder, string graph_file): folder(folder), graph_file(graph_file)
    {
		srand(time(NULL));
		sfmt_init_gen_rand(&sfmtSeed, rand());

        readNM();

		for (int i = 0; i < n*2; i++)
        {

			OutNode.push_back(vector<int>());
            InNode.push_back(vector<int>());
			hasnode.push_back(false);
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

