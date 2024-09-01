//
// Created by DBL-XQ on 2024/3/19.
//

//#define HEAD_TRACE
#define HEAD_INFO

#define HEAD_INFO
//#define HEAD_TRACE
#include "sfmt/SFMT.h"
#include "head.h"
#include "memoryusage.h"
#include "graph_tim.h"

void run(TimGraph & m, string dataset, string seedset, int k, double epsilon, string model, int p ){
    cout << "dataset:" << dataset << " k:" << k << " p:" << p << " epsilon:"<< epsilon <<   " model:" << model << endl;
    m.k=k;
    if(model=="IC")
        m.setInfuModel(InfGraph::IC);
    else if(model=="LT")
        m.setInfuModel(InfGraph::LT);
    else
    ASSERT(false);

    cout<<"Finish Read Graph, Start Influecne Maximization"<<endl;
    double tau=m.EstimateOPT(epsilon);
    cout<<"Time used: " << Timer::timeUsed[100]/TIMES_PER_SEC << "s" <<endl;
    cout<<"tau="<<tau<<endl;
    cout<<"Selected k SeedSet: ";
    for(auto item:m.seedSet)
        cout<< item << " ";
    cout<<endl;
    cout<<"Estimated Influence: " << m.InfluenceHyperGraph() << endl;
    Counter::show();

    ////ouput
    if (p ==100){
        string snum = std::to_string(k);
        string outfile = seedset +"seed_IM_"+snum+".txt";
        ofstream of(outfile);
        for(auto item:m.seedSet)
            of<< item << endl;
        of.close();
    }
    else{
        string per = "_p"+to_string(p);
        string snum = std::to_string(k);
        string outfile = seedset +"seed_IM_"+snum+per+".txt";
        cout << outfile <<endl;
        ofstream of(outfile);
        for(auto item:m.seedSet)
            of<< item << endl;
        of.close();
    }

}
void parseArg(int argn, char ** argv)
{
    string dataset="";
    string seedset="";

    double epsilon=0;
    string model="";
    int k=0;
    int p=100;

    for(int i=0; i<argn; i++)
    {// -dataset ../dataset/Email-Core/ -model IC  -seedset ../seedset/Email-Core/  -epsilon 0.1 -k 50
        if(argv[i]==string("-dataset")) dataset=string(argv[i+1]);
        if(argv[i]==string("-seedset")) seedset=string(argv[i+1]);
        if(argv[i]==string("-epsilon")) epsilon=atof(argv[i+1]);
        if(argv[i]==string("-k")) k=atoi(argv[i+1]);
        if(argv[i]==string("-p")) p=atoi(argv[i+1]);
        if(argv[i]==string("-model")) {
            if(argv[i+1]==string("LT"))
            {
                model=argv[i+1];
            }
            else if(argv[i+1]==string("IC"))
            {
                model=argv[i+1];
            }
            else
                ExitMessage("model should be IC or LT");
        }
    }
    if (dataset=="")
        ExitMessage("argument dataset missing");
    if (k==0)
        ExitMessage("argument k missing");
    if (epsilon==0)
        ExitMessage("argument epsilon missing");
    if (model=="")
        ExitMessage("argument model missing");


    string graph_file;
    string per = to_string(p);
    if (p == 100){
        graph_file=dataset + "graph.txt";
    }
    else{
        graph_file=dataset + "graph_p"+per+".txt";
    }
//    if(model=="IC")
//        graph_file=dataset + "graph_ic.inf";
//    else if(model=="LT")
//        graph_file=dataset + "graph_lt.inf";

    TimGraph m(dataset, graph_file, per);
    run(m, dataset, seedset, k ,  epsilon, model, p );
    cout<<"second time"<<endl;
    disp_mem_usage("");
    cout<<endl;
}





int main(int argn, char ** argv)
{
    OutputInfo info(argn, argv);
    parseArg( argn, argv );
}