//
// Created by DBL-XQ on 2024/3/26.


#include "InfCompute.h"
#include <omp.h>

set<int> readseedset(string seedfile, const int seednum);

int main(int argn, char*argv[])
{
    if (argn == 1){//
        cout << "Usage: ./InfExpectation_computation --theta 10000 --dataset ../dataset/Pokec/ --model IC --mode IM --seedset ../seedset/Pokec/  --seednum 50" << endl;
        exit(0);
    }

    string model;
    string dataset;
    string seedset;
    string result;
    string mode; // random IM inf
    int theta = 10000;
    int seednum = 10;
    int p = 100;

    for (int i = 1; i<argn; i++){
        if (string(argv[i]) == "--theta")
            theta = atoi(argv[i + 1]);
        if (string(argv[i]) == "--seednum")
            seednum = atoi(argv[i + 1]);
        if (string(argv[i]) == "--p")
            p = atoi(argv[i + 1]);
        if (string(argv[i]) == "--seedset")
            seedset = argv[i + 1];
        if (string(argv[i]) == "--dataset")
            dataset = argv[i + 1];
        if (string(argv[i]) == "--result")
            result = argv[i + 1];
        if (string(argv[i]) == "--model")
            model = argv[i + 1];
        if (string(argv[i]) == "--mode")
            mode = argv[i + 1];
    }

    string diffusionPath;
    string per = to_string(p);
    string seedfile;
    string snum = std::to_string(seednum);
    if (p == 100){
        diffusionPath = dataset + "graph.txt";
        seedfile = seedset + "seed_" + mode +"_" + snum + ".txt";
    }
    else{
        diffusionPath=dataset + "graph_p"+per+".txt";
        seedfile = seedset + "seed_" + mode +"_" + snum + "_p" + per + ".txt";
    }

    InfC mi(dataset, diffusionPath, per);

    if (model == "IC"){
        mi.setInfuModel(Graph::IC);
    }
    else if (model == "IC_van"){
        mi.setInfuModel(Graph::IC_van);
    }
    else if (model == "LT"){
        mi.setInfuModel(Graph::LT);
    }


    mi.seedSet = readseedset(seedfile,seednum);

        ////////////////////////////////// compute influenced expectation /////////////////////////////////
        double startTime = omp_get_wtime();
        vector<double> influence = mi.Influence(theta);
        double endTime = omp_get_wtime();

//    set <int> :: iterator iter1;
//    for(iter1 = mi.nodes.begin(); iter1 != mi.nodes.end(); iter1 ++){
//        auto nodeId = *iter1;
//        cout << "node id " << nodeId << "\t influence expectation: " << influence[nodeId] << endl;
//    }


        cout << "Number of seeds: " << mi.seedSet.size() << endl;
        cout << "theta: " << theta << endl;
        //cout << "startTime: " << startTime << "\t endTime: " << endTime<< endl;
        cout << "Total Time: " << (double) (endTime - startTime) << " s" << endl;




        //// write in outfile the influenced expectation Ie(v,S) of each node in dataset
        if (p == 100){
            string outfile;
            outfile = dataset + "inf_exp_" + model +"_" + mode + "_s" + snum + ".txt";
            ofstream out(outfile);
            if (!out.is_open())
                std::cout << "The file \"" + outfile + "\" can NOT be opened\n";
            set<int>::iterator iter;
            for (iter = mi.nodes.begin(); iter != mi.nodes.end(); iter++) {
                auto nodeId = *iter;
                //cout << "node id " << nodeId << "\t influence expectation: " << influence[nodeId] << endl;
                out << nodeId << " " << influence[nodeId] << endl;
            }
            out.close();

            string timefile;
            //timefile = dataset + "time_inf_exp_" + model +"_" + mode + "_s" + snum + ".txt";
            timefile = dataset + "time_inf_exp_IC_" + mode + "_s" + snum + ".txt";
            ofstream timeout(timefile);
            if (!timeout.is_open())
                std::cout << "The file \"" + timefile + "\" can NOT be opened\n";
            timeout << "theta: " << theta << endl;
            timeout << "model: " << model << endl;
            timeout << "seed_num: " << seednum << endl;
            timeout << "Total influence expectation computation time (s): " << (double)(endTime - startTime)  << endl;
            timeout.close();
        }
        else{
            string outfile;
            per = "_p"+to_string(p);
            outfile = result + "inf_exp_" + model +"_" + mode + "_s" + snum + per + ".txt";
            ofstream out(outfile);
            if (!out.is_open())
                std::cout << "The file \"" + outfile + "\" can NOT be opened\n";
            set<int>::iterator iter;
            for (iter = mi.nodes.begin(); iter != mi.nodes.end(); iter++) {
                auto nodeId = *iter;
                //cout << "node id " << nodeId << "\t influence expectation: " << influence[nodeId] << endl;
                out << nodeId << " " << influence[nodeId] << endl;
            }
            out.close();
        }


}


set<int> readseedset(string seedfile, const int seednum)
{
    cout << "Read seedfile ..." << endl;
    std::ifstream infile(seedfile);
    if (!infile.is_open()) {
        cout<< "Cannot open file " << seedfile << endl;
    }
    set<int> s;
    int seed;
    cout <<"seedset["<< seednum<< "]: ";
    while (infile >> seed) {
        s.insert(seed);
        cout << seed <<"\t";
    }
    cout <<endl;
    infile.close();
    return s;
}

