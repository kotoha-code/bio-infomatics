#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <random>
using namespace std;

#define NUM_FEATURES
#define NUM_SEQS

vector<string>LoadSolubilityFile(string filename,vector<string> feature_name, vector<vector<double>> dataset, vector<int> labels){
    ifstream fin(filename);
    vector<string> seqs;
    string line;
}

vector<string>DivideDataset(vector<vector<double>> dataset,vector<int> labels, vector<vector<double>> training_dataset,vector<int> training_labels,vector<vector<double>> test_dataset,vector<int> test_labels,double test_ratio){

}

int main(void){
    vector<string> feature_name(NUM_FEATURES, "");
    vector<vector<double>> dataset(NUM_SEQS, vector<double>(NUM_FEATURES,0.0));
    vector<int>labels(NUM_SEQS);

    LoadSolubilityFile("protein_solubility_dataset.txt",feature_name,dataset,labels);

    vector<vector<double>>training_dataset;
    vector<int>training_labels;
    vector<vector<double>>test_dataset;
    vector<int>test_labels;
    double test_ratio=0.2;

    DivideDataset(dataset,labels,training_dataset,training_labels,test_dataset,test_labels,test_ratio);
}