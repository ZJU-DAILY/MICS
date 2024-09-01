//
// Created by DBL-XQ on 2024/3/5.
//

#ifndef INC_20240217_INFEXPECTATION_COMPUTATION_OPTION_H
#define INC_20240217_INFEXPECTATION_COMPUTATION_OPTION_H

#endif //INC_20240217_INFEXPECTATION_COMPUTATION_OPTION_H

#ifndef _OPTION_H_
#define _OPTION_H_
#include <sstream>
#include <map>

class OptionParser{

private:
    std::map<std::string, char *> paras;
    bool isValid;
public:
    OptionParser(int argc, char ** argv);
    char * getPara(const char * str);
    bool validCheck();
};

#endif
