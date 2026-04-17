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
    //ファイルの取得、名前と配列別々に
    vector<string> names;
    vector<string> sequences;
    string line;

    ifstream fin("data/seq/REB1.promoters");
    while (getline(fin, line)) {
        if (line[0] == '>') {
            names.push_back(line.substr(1));
            sequences.push_back(""); 
        } else {
            sequences.back() += line; 
        }
    }
    map<char,int> idx = {{'A',0},{'C',1},{'G',2},{'T',3}};

    //1.各配列に対して、ランダムにモチーフの出現位置を決定して、それを初期値とする。
    mt19937 rng(42);
    vector<int> instance; //初期値を格納する配列
    for(int i=0;i<sequences.size();i++){
        uniform_int_distribution<int> dist(0, sequences[i].size()); 
        int r = dist(rng);
        instance[i]=r;
    }

    for(int i=0;i<10000000;i++){
    //2.N本の配列群から配列を1つランダムに選び、残りのN−1本の配列群に対して、モチーフの出現位置からスコア行列を算出する。

    }



}