#include "sbml.hpp"

VectorXd simulateSBML(int useDet, double ti, double tf, const VectorXd &c0, const VectorXd &k){
    RoadRunner r = RoadRunner("sbml/model_sbml.xml");
    vector<double> init = convertInit(c0);
    r.changeInitialConditions(init);
    SimulateOptions opt;
    opt.start = ti;
    opt.duration = tf;
    opt.steps = 2;

    // to apparently change the values of parameters in the model, we must first feed the vector into a double array.
    double kVals[k.size()];
    for(int i = 0; i < k.size(); ++i){
        kVals[i] = k(i);
    }
    r.getModel()->setGlobalParameterValues(k.size(),0,kVals); // set new global parameter values here.
    if(useDet > 0){
        r.setIntegrator("cvode");
    }else{
        r.setIntegrator("gillespie");
    }
    const DoubleMatrix res = *r.simulate(&opt);
    VectorXd evolved = VectorXd::Zero(res.numCols() - 1);
    for(int i = 1; i < res.numCols(); i++){
        evolved(i - 1) = res[res.numRows() - 1][i];
    }
    return evolved;
}

vector<string> getSpeciesNames(const string& path){
    vector<string> listOfSpecies;
    tinyxml2::XMLDocument doc;
    doc.LoadFile(path.c_str());
    tinyxml2::XMLElement* xmlModel = doc.FirstChildElement() -> FirstChildElement("model");
    for (tinyxml2::XMLElement* child = xmlModel->FirstChildElement(); child != NULL; child = child->NextSiblingElement())
    {
        if(strcmp(child->Value(), "listOfSpecies") == 0){
            for (tinyxml2::XMLElement* species = child->FirstChildElement(); species != NULL; species = species->NextSiblingElement())
            {
                const char* name;
                species->QueryStringAttribute("name", &name);
                listOfSpecies.push_back(std::string(name));
            }
        }
    }
    return listOfSpecies;
}
vector<int> specifySpeciesFromProteinsList(const string& path, vector<string> &species, int nObs){
    std::ifstream input(path);
    vector<string> pro;
    vector<int> indices;
    if(!input.is_open()){
        throw std::runtime_error("Could not open proteins file!");
        exit;
    }
    string line;
    while(std::getline(input,line)){
        pro.push_back(line);
    }
    for(int i = 0; i < pro.size(); i++){
        string poi = pro[i];
        for(int j = 0; j < species.size(); j++){
            if(poi == species[j]){
                indices.push_back(j);
            }
        }
    }
    if(nObs != indices.size()){
        cout << "Error Mismatch in Number of Columns (Species) Provided in Data Files Versus Proteins Specified in " << path << endl;
        cout << "Observed:" << nObs << " Number in POI:" << indices.size() << endl; 
        exit(EXIT_FAILURE);
    }
    input.close();
    return indices;
}
