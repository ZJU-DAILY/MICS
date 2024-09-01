#include "stdafx.h"

int main(int argc, char * argv[]) {


    if (argc == 1){
        cout << "Usage: ./InfCS -dataset=../../dataset/test/ -model=IC -mode=random -func=avg "
                "-alg=approx -med=global_naive -seed=2 -k=1 -l=1 -s=2" << endl;
        exit(0);
    }

    const TArgument Arg(argc, argv);

    SMIC *smic = new SMIC();

    //// ======================================= read graph ============================================
    string attribute_file, graph_file, inf_file;
    attribute_file = Arg._dataset + "attribute.txt";
    graph_file = Arg._dataset + "graph.txt";
    string snum = std::to_string(Arg._seednum);
    if (Arg._inf == "SFIS"){
        inf_file = Arg._dataset + "inf_exp_" + Arg._model  + "_"+ Arg.seed_mode + "_s" + snum + ".txt";
    }
    else if (Arg._inf == "ONPR"){
        inf_file = Arg._dataset + "preprob_RF_" + Arg.seed_mode + "_" + snum + ".txt";
    }
    if (Arg._p == 100){
        smic->readNM(attribute_file);
        smic->readGraph(graph_file);
        smic->readInf(inf_file);
    }
    else{
        // test scalability
        string p = std::to_string(Arg._p);
        attribute_file = Arg._subgraph + "attribute_p" + p +".txt";
        graph_file = Arg._subgraph + "graph_p" + p +".txt";
        if (Arg._inf == "SFIS"){
            inf_file = Arg._dataset + "inf_exp_" + Arg._model  + "_"+ Arg.seed_mode + "_s" + snum +"_p"+p + ".txt";
        }
        else if (Arg._inf == "ONPR"){
            inf_file = Arg._dataset + "preprob_RF_" + Arg.seed_mode + "_" + snum +"_p"+p + ".txt";
        }
        smic->readNM(attribute_file);
        smic->readsubGraph(graph_file, Arg._p);
        smic->readsubInf(inf_file, Arg._p);
    }

    smic->graphinf();

    //// ======================================= Algorithm ============================================
    if (Arg._func == "avg") {
        if (Arg._algName == "exact") {
            if (Arg._method == "global") {
                //smic->basic_random_avg_global_mic(Arg._dataset, Arg._k, Arg._l,Arg._s);
                smic->basic_minweight_avg_global_mic(Arg._dataset, Arg._k, Arg._l,Arg._s,Arg._inf,snum);
            }
            else if (Arg._method == "local") {
                smic->basic_maxweight_avg_local_mic(Arg._dataset, Arg._k, Arg._l,Arg._s,Arg._inf,snum);
            }
        }
        else if (Arg._algName == "approx") {
            if (Arg._method == "global") { /////////////////////////////////////////// global
                smic->improved_minweight_avg_global_mic(Arg._dataset, Arg._k, Arg._l,Arg._s,Arg._inf,snum, Arg._p);
                }
            else if (Arg._method == "local") { /////////////////////////////////////////// local
                smic->improved_maxweight_avg_local_mic(Arg._dataset, Arg._k, Arg._l,Arg._s,Arg._inf,snum, Arg._p);
            }
        }
    }
    else if (Arg._func == "min") {
        if (Arg._method == "global") {
            smic->improved_minweight_min_global_mic(Arg._dataset, Arg._k, Arg._l,Arg._s,Arg._inf,snum);
        }
        else if (Arg._method == "local") {
            smic->improved_one_maxweight_min_local_mic(Arg._dataset, Arg._k, Arg._l,Arg._s,Arg._inf,snum);
        }
    }
    else if (Arg._func == "max") {
        if (Arg._method == "global") {
            smic->improved_minweight_max_global_mic(Arg._dataset, Arg._k, Arg._l,Arg._s,Arg._inf,snum);
        }
        else if (Arg._method == "local") {
            smic->improved_maxweight_max_local_mic(Arg._dataset, Arg._k, Arg._l,Arg._s,Arg._inf,snum);
        }
    }
    else if (Arg._func == "sum") {
        smic->improved_sum_mic(Arg._dataset, Arg._k, Arg._l,Arg._s,Arg._inf,snum);
    }


    std::cout << "\n============================ End ==============================" << std::endl;
    return 0;
}
