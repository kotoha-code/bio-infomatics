#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <random>
using namespace std;

//一行ずつファイルを読み込んでseqsというファイルにそのまま写しとる関数
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
  
//MATa1を読み込む
    ifstream ist("data");
    if(!ist){
       cerr << "Cannot open data" << endl;
       exit(1);
    }
//塩基と対応するMapを作る
    map<char,int> idx = {{'A',0},{'C',1},{'G',2},{'T',3}};

   vector<string> filename= {"MATa1", "MATalpha2", "MCM1", "MIG1", "PHO4", "RCS1", "ROX1", "TAF"};

   for(int name=0;name<8;name++){
        vector<string> seqs = read_motif("data/motif/"+filename[name]);
//filenameを紹介
    cout<<"filename : "<<filename[name]<<endl;
    //ファイル一行ずつに含まれる文字数
    int L = seqs[0].size();
    
    //頻度表を作成
    vector<vector<double>> count(4, vector<double>(L,1.0)); 

    for(int i=0;i<seqs.size();i++){
        for(int j=0;j<L;j++){
            count[idx[seqs[i][j]]][j]++;
        }
    }

/*//一旦出力
    for(int i=0;i<seqs.size();i++){
        for(int j=0;j<L;j++){
            cout<<count[i][j]<<" ";
        }
        cout<<endl;
    }*/

// pを計算
    vector<vector<double>> p(4, vector<double>(L));
    double bunbo=seqs.size()+4;
     for(int i=0;i<p.size();i++){
        for(int j=0;j<L;j++){
            p[i][j]=count[i][j]/bunbo;
        }
    }
    /*//一旦出力
     for(int i=0;i<seqs.size();i++){
        for(int j=0;j<L;j++){
            cout<<p[i][j]<<" ";
        }
        cout<<endl;
    }*/

// q(x)を計算した
    double sum=7519429+7519429+4637676+4637676;
    vector<double> q{7519429/sum,4637676/sum,4637676/sum,7519429/sum};
//オッズスコア計算    
   vector<vector<double>> score(4, vector<double>(L));
     for(int i=0;i<score.size();i++){
        for(int j=0;j<L;j++){
            score[i][j]=log(p[i][j]/q[i]);
        }
    }
    //一旦出力
     for(int i=0;i<score.size();i++){
        for(int j=0;j<L;j++){
            cout<<score[i][j]<<" ";
        }
        cout<<endl;
    }

//課題２
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
        vector<vector<double>> hitv(sequences.size(), vector<double>(sequences[0].size()-L,0.0)); 

        for(int num=0;num<sequences.size();num++){
           for(int i=0; i<sequences[num].size()-L;i++){
                double hit=0.0;
                for(int j=0;j<L;j++){
                    hit+=score[idx[sequences[num][i+j]]][j];
                }
                hitv[num][i]=hit;
           }
        }
        for(int i=0;i<hitv.size();i++){
            cout<<names[i]<<" "<<endl;
            for(int j=0;j<hitv[0].size();j++){
                if(hitv[i][j]>5.0){
                    cout<<"pos : "<< i+j<<","<<hitv[i][j]<<" "<<endl;
                }
            }
            cout<<endl;
        }

        //課題３: 閾値の決定
    // バックグラウンド確率q(x)に従うランダム配列を大量生成し
    // スコアの分布からp値に対応する閾値を求める

    const int N_RANDOM = 10000;  // ランダム配列の本数
    const int SEQ_LEN  = 500;    // プロモーターと同じ長さ
    const double P_VALUE = 0.001; // p値

    // ランダム配列生成の準備
    // q[0]=A, q[1]=C, q[2]=G, q[3]=T の確率で塩基を選ぶ
    mt19937 rng(42);
    discrete_distribution<int> base_dist({q[0], q[1], q[2], q[3]});
    const string BASES = "ACGT";

    // 全ランダム配列のスコアを集める
    vector<double> all_scores;
    all_scores.reserve((long long)N_RANDOM * (SEQ_LEN - L + 1));

    for (int r = 0; r < N_RANDOM; r++) {
        // ランダム配列を1本生成
        string rand_seq(SEQ_LEN, 'A');
        for (int i = 0; i < SEQ_LEN; i++)
            rand_seq[i] = BASES[base_dist(rng)];

        // その配列の全位置でスコアを計算
        for (int i = 0; i <= SEQ_LEN - L; i++) {
            double hit = 0.0;
            for (int j = 0; j < L; j++)
                hit += score[idx[rand_seq[i+j]]][j];
            all_scores.push_back(hit);
        }
    }

    // スコアを昇順にソート
    sort(all_scores.begin(), all_scores.end());

    // 上位p値の割合に対応するスコアが閾値
    int threshold_idx = (int)(all_scores.size() * (1.0 - P_VALUE));
    double threshold = all_scores[threshold_idx];

    cout << "閾値 (p=" << P_VALUE << "): " << threshold << endl;

    // 閾値を使って結合部位を判定
    cout << "=== 結合部位予測 (閾値=" << threshold << ") ===" << endl;
    for (int num = 0; num < sequences.size(); num++) {
        for (int i = 0; i < hitv[num].size(); i++) {
            if (hitv[num][i] >= threshold) {
                cout << filename[name] << " -> " << names[num]
                     << "  pos=" << i
                     << "  score=" << hitv[num][i] << endl;
            }
        }
    }
  } 
  return 0; 
}
