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
    cout << "Please Make Sure to Specify \"-x <DIRECTORY>\" in call, otherwise it will default to data/X directory, which may be inside of the docker container and not on your device!" << endl;
    return "data/X";
}
string getYPath(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-y");
    if(flag != -1){
        return argv[flag+1];
    }
    cout << "Please Make Sure to Specify \"-y <DIRECTORY>\" in call, otherwise it will default to data/Y directory, which may be inside of the docker container and not on your device!" << endl;
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
        cout << "To specify model path, do: ./BNGMM -m <path> i.e ./BNGMM -m model.bngl" << endl
        << "To specify config path, do: ./BNGMM -c <path> i.e ./BNGMM -c Config.csv" << endl
        << "To specify time steps file path, do: ./BNGMM -t <path> i.e ./BNGMM -t time_steps.csv" << endl
        << "To specify simulation truth rate constants path, do: ./BNGMM -r <path> i.e ./BNGMM -r true_rates.csv" << endl
        << "To specify data X directory where X data files are located, do: ./BNGMM -x <DIRECTORY> i.e ./BNGMM -x to/X/dir" << endl
        << "To specify data Y directory where true Y data files are located, do: ./BNGMM -y <DIRECTORY> i.e ./BNGMM -y to/Y/DIRECTORY" << endl
        << "To specify an output directory where output files such as graphing and output txt files, ./BNGMM -o <path> i.e ./BNGMM -o /frontend/graphs/6pro" << endl
        << "If you have more species in the system than observed protein species, then please supply a list of proteins in a .txt file." << endl
        << "i.e ./BNGMM -p listOfObservedProteinsInOrder.txt " << endl;
        return true;
    }
    return false;
}

bool holdRates(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-hr");
    return flag != -1;
}

bool seedRates(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-s");
    return flag != -1;
}

bool forecast(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-f");
    return flag != -1;
}

string getHeldRatesDir(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-hr");
    return argv[flag+1];
}

string getSeededRates(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-s");
    return argv[flag+1];
}

string getForecastedTimes(int argc, char **argv){
    int flag = getIndexFlag(argc, argv, "-f");
    return argv[flag+1];
}
