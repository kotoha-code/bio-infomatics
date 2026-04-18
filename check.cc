#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <random>
using namespace std;

int main(void){
 ifstream ist("data");
    if(!ist){
       cerr << "Cannot open data" << endl;
       exit(1);
    }
//塩基と対応するMapを作る
    map<char,int> idx = {{'A',0},{'C',1},{'G',2},{'T',3}};

    vector<double> count(4,0);

    vector<string> names;
    vector<string> sequences;
    string line;
    ifstream fin("data/seq/promoters");
    while (getline(fin, line)) {
        if (line[0] == '>') {
            names.push_back(line.substr(1));
            sequences.push_back(""); 
        } else {
            sequences.back() += line; 
        }
    }
    int L = sequences[0].size();

    for(int i=0;i<sequences.size();i++){
        for(int j=0;j<L;j++){
            count[idx[sequences[i][j]]]++;
        }   
    }    
    double total;
    for(int i=0;i<count.size();i++){
        cout<<count[i]<<" ";
        total+=count[i];
    }
    cout<<endl;
    for(int i=0;i<count.size();i++){
        //バックグラウンド出現頻度
        cout<<count[i]/total<<" ";
    }
    cout<< endl;
    return 0;
}