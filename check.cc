#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <random>
using namespace std;


vector<string> read_motif(string filename){
    ifstream fin(filename);
    vector<string> seqs;
    string line;

    while(getline(fin,line)){
        if(line.size()==0) continue;
        seqs.push_back(line);
    }
    return seqs;
}

int main(void){
 ifstream ist("data");
    if(!ist){
       cerr << "Cannot open data" << endl;
       exit(1);
    }
//塩基と対応するMapを作る
    map<char,int> idx = {{'A',0},{'C',1},{'G',2},{'T',3}};

    vector<string> filename= {"MATa1", "MATalpha2", "MCM1", "MIG1", "PHO4", "RCS1", "ROX1", "TAF"};
    vector<double> count(4,0);

    for(int name=0;name<8;name++){
        vector<string> seqs= read_motif("data/motif/"+filename[name]);
        int L = seqs[0].size();

        for(int i=0;i<seqs.size();i++){
            for(int j=0;j<L;j++){
                count[idx[seqs[i][j]]]++;
            }   
        }
    }
    double total;
    for(int i=0;i<count.size();i++){
        //全モチーフに含まれている塩基数（ATGC順）
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