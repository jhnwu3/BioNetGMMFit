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
string getOutputPath(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-o");
    if(flag != -1){
        return argv[flag+1];
    }
    return "/frontend/graph";
}
bool outPathExists(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-o");
    return flag != -1;
}
bool modelPathExists(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-m");
    return flag != -1;
}

bool graphingEnabled(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "--g");
    return flag != -1;
}

string getProPath(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-p");
    return argv[flag+1];
}
bool proPathExists(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-p");
    return flag != -1;
}
bool helpCall(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-h");
    if(flag != -1){
        cout << "To specify model path, do: ./CyGMM -m <path> i.e ./CyGMM -m model.bngl" << endl
        << "To specify config path, do: ./CyGMM -c <path> i.e ./CyGMM -c Config.csv" << endl
        << "To specify time steps file path, do: ./CyGMM -t <path> i.e ./CyGMM -t time_steps.csv" << endl
        << "To specify simulation truth rate constants path, do: ./CyGMM -r <path> i.e ./CyGMM -r true_rates.csv" << endl
        << "To specify data X directory where X data files are located, do: ./CyGMM -x <path> i.e ./CyGMM -x model.bngl" << endl
        << "To specify data Y directory where true Y data files are located, do: ./CyGMM -y <path> i.e ./CyGMM -y model.bngl" << endl
        << "To specify an output directory where output files such as graphing and output txt files, ./CyGMM -o <path> i.e ./CyGMM -o /frontend/graphs/6pro" << endl
        << "If you have more species in the system than observed protein species, then please supply a list of proteins in a .txt file." << endl
        << "i.e ./CyGMM -p listOfObservedProteinsInOrder.txt " << endl;
        return true;
    }
    return false;
}