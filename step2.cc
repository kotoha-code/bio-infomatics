#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
using namespace std;

//一行ずつファイルを読み込んでseqsというファイルにそのまま写しとる関数
vector<string> read_motif(string filename){
    ifstream fin(filename);
    vector<string> seqs;
    vector<string> seqs_only_base;
    string line;

    while(getline(fin,line)){
        if(line.size()==0) continue;
        seqs.push_back(line);
    }
    return seqs;
}