//
// Created by DBL-XQ on 2024/3/5.
//

#include "option.h"
#include "graph.h"


#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <fstream>
#include <cstring>
#include <random>

//// -n 1005 -o ../../datasets/EmailCore/fake.seeds -k 10 -m top -f 0.02 -s ../../datasets/EmailCore/inf.txt

using namespace std;

float estimateFakeInfluenceParallel(Graph &g, int n);


void generateSeedsTop(int n, int k, double f, string infFile, vector<int> &seeds) {
    ifstream in(infFile);
    int i, v1, rand_pos;
    float v2;
    int num = n * f;
    vector<int> candidates(num,0);
    srand (time(NULL));

    for (i = 0; i < num; i++) {
        in >> v1 >> v2;
        rand_pos = rand() % (i+1);
        if (rand_pos != i) candidates[i] = candidates[rand_pos];
        candidates[rand_pos] = v1;
    }
    in.close();

    for (i = 0; i < k; i++) {
        seeds[i] = candidates[i];
    }
}

void generateSeedsRandom(int n, int k, vector<int> &seeds) {
    int i, v, rand_pos;
    vector<int> candidates(n,0);
    srand (time(NULL));

    for (i = 0; i < n; i++) {
        candidates[i] = i;
    }

    for (i = n-1; i > 0; i--) {
        rand_pos = rand() % (i+1);
        v = candidates[rand_pos];
        candidates[rand_pos] = candidates[i];
        candidates[i] = v;
    }

    for (i = 0; i < k; i++) {
        seeds[i] = candidates[i];
    }
}

void generateSeedsComm(Graph &g, int n, int k, const char* commFile, vector<int> &seeds) {
    return;
}

int main(int argc, char ** argv)
{

    srand(time(NULL));

    if (argc == 1){ // --method random --k 100 --f 0.2 --dataset ../dataset/Pokec/ --seedset ../seedset/Pokec/ --r 1000 --model IC
        cout << "Usage: ./InfExpectation_computation --method inf --k 20 --f 40 --dataset ../dataset/Wiki-Vote/ --seedset ../seedset/Wiki-Vote/" << endl;
        exit(0);
    }

    string seedset; // result seedset path
    string dataset; // node_inf.txt
    int k; // number of seeds
    double f; // method inf top-ratio
    const char *method; // random/inf


    for (int i = 1; i<argc; i++){
        if (string(argv[i]) == "--k")
            k = atoi(argv[i + 1]);
        if (string(argv[i]) == "--f")
            f = atof(argv[i + 1]);
        if (string(argv[i]) == "--dataset")
            dataset = argv[i + 1];
        if (string(argv[i]) == "--seedset")
            seedset = argv[i + 1];
        if (string(argv[i]) == "--method")
            method = argv[i + 1];
    }



    string infFile = dataset + "node_inf_IC.txt";

    vector<int> seeds(k,0);

    string dir,filename;
    dir = dataset;
    filename = dataset + "graph.txt";
    Graph g(dir, filename);
    int n = g.n;

    if (strcmp(method, "inf") == 0) {
        cout << "generateSeedsTop ..." <<endl;
        generateSeedsTop(n, k, f, infFile, seeds);
    } else if (strcmp(method, "random") == 0) {
        cout << "generateSeedsRandom ..." <<endl;
        generateSeedsRandom(n, k, seeds);
    } else {
        printf("Incorrect method option!");
        return -1;
    }

    string snum = std::to_string(k);
    string outfile = seedset +"seed_"+method+"_"+snum+".txt";
    ofstream of(outfile);
    //of << k << endl;
    for (int i = 0; i < k; i++){
        of << seeds[i] << endl;
        cout << seeds[i] << " ";
    }cout << endl;
    of.close();

    //
    // ALL DONE
    //
}