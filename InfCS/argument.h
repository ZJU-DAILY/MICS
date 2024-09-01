#ifndef WEIGHT_COMMUNITY_ARGUMENT_H
#define WEIGHT_COMMUNITY_ARGUMENT_H

class Argument
{
public:
    int _k = 5; // out_degree constraint k.
    int _l = 5; // in_degree constraint l.
    int _s = 20; // size constraint.
    int _seednum = 10; // size of seedset
    int _p = 100; // percentage of nodes in graph
    std::string _subgraph = "../../subgraphset/Pokec/";
    std::string _dataset = "../../dataset/Wiki-Vote/";
    std::string _model = "IC"; // propagation model. IC, LT
    std::string seed_mode = "IM"; // seed_mode IM random inf
    std::string _algName = "approx"; // Algorithm. Default is approx. exact, approx
    std::string _func = "avg"; // Aggregated function. avg, min, max, sum
    std::string _method = "global"; // Method. global, local
    std::string _inf = "SFIS"; // Inf Comp method. SFIS, ONPR

    Argument(int argc, char * argv[])
    {
        std::string param, value;
        for (int ind = 1; ind < argc; ind++)
        {
            if (argv[ind][0] != '-') {
                break;
            }
            std::stringstream sstr(argv[ind]);
            getline(sstr, param, '=');
            getline(sstr, value, '=');
            if (!param.compare("-k")) {
                _k = stoi(value);
            }
            if (!param.compare("-l")) {
                _l = stoi(value);
            }
            else if (!param.compare("-s")) {
                _s = stoi(value);
            }
            else if (!param.compare("-seed")) {
                _seednum = stoi(value);
            }
            else if (!param.compare("-p")) {
                _p = stoi(value);
            }
            else if (!param.compare("-subgraph")) {
                _subgraph = value;
            }
            else if (!param.compare("-dataset")) {
                _dataset = value;
            }
            else if (!param.compare("-model")) {
                _model = value;
            }
            else if (!param.compare("-mode")) {
                seed_mode = value;
            }
            else if (!param.compare("-alg")) {
                _algName = value;
            }
            else if (!param.compare("-func")) {
                _func = value;
            }
            else if (!param.compare("-med")) {
                _method = value;
            }
            else if (!param.compare("-inf")) {
                _inf = value;
            }
        }
    }
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;

#endif //WEIGHT_COMMUNITY_ARGUMENT_H
