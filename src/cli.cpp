#include "cli.hpp"

int getIndexFlag(int argc, char** argv, const string &flag){
    for(int i = 0; i < argc; ++i){
        if(argv[i] == flag){ // return index of flag
            return i;
        }
    }
    return -1; // return - 1 if can't find specified flag
}

string getConfigPath(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-c");
    if(flag != -1){
        return argv[flag+1];
    }
    return "Config.csv";
}
string getTimeStepsPath(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-t");
    if(flag != -1){
        return argv[flag+1];
    }
    return "time_steps.csv";
}
string getTrueRatesPath(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-r");
    if(flag != -1){
        return argv[flag+1];
    }
    return "true_rates.csv";
}
string getXPath(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-x");
    if(flag != -1){
        return argv[flag+1];
    }
    return "data/X";
}
string getYPath(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-y");
    if(flag != -1){
        return argv[flag+1];
    }
    return "data/Y";
}
string getModelPath(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-m");
    if(flag != -1){
        return argv[flag+1];
    }
    return "model.bngl";
}
