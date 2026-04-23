#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
using namespace std;

#define NUM_FEATURES 53
#define NUM_SEQS 10000

struct TreeNode {
    int    feature_id;    // 分岐に使う特徴量のインデックス
    double threshold;     // 分岐閾値（以下なら左、超えたら右）
    int    left_class_id; // 左に落ちたデータの予測ラベル
    int    right_class_id;// 右に落ちたデータの予測ラベル
};

void LoadSolubilityFile(string filename,vector<string>& feature_name, vector<vector<double>>& dataset, vector<int>& labels){
    ifstream ifs(filename);
    if(!ifs){
        cerr << "Cannot open data" << endl;
        return;
    }
    //ifsから1行読んでlineに入れる
    string line;
    getline(ifs,line);
    istringstream ss(line); //lineをストリームとして扱う文法らしい
    string temp; //取ってきた一単語をとりあえず保存するために
    vector<string> headers; 

    //データの1行目を区切ってheaderに入れていく
    while(getline(ss,temp,' ')){
        headers.push_back(temp);
    }
    cout<<headers.size()<<endl;
    //headerの中でid名以外を格納
    for(int i=1; i<=NUM_FEATURES; i++){
        feature_name[i-1]=headers[i];
    }

    int row=0;
    while(getline(ifs,line)&&row<NUM_SEQS){
        istringstream ss(line); //読んだ行をstreamに変換
        getline(ss,temp,',');  //protein idは使わないから捨てる

        for(int i=0;i<NUM_FEATURES;i++){
            getline(ss,temp,',');
            dataset[row][i]=stod(temp); //stod=文字列から少数に変換する
        }
        getline(ss,temp,',');
        labels[row]=stoi(temp); //stoi=文字列から整数に変換する
        row++;
    }
}

vector<string>DivideDataset(vector<vector<double>>& dataset,vector<int>& labels, vector<vector<double>>& training_dataset,vector<int>& training_labels,vector<vector<double>>& test_dataset,vector<int>& test_labels,double test_ratio){
    //インデックス配列作る
    vector<int>indices(dataset.size());
    iota(indices.begin(),indices.end(),0);
    //シャッフル
    mt19937 gen(0);
    shuffle(indices.begin(),indices.end(),gen);

    //データセットの先頭20%をテスト用に格納
    int testsize=(int)(dataset.size()*test_ratio);
    for(int i=0;i<testsize;i++){
        test_dataset.push_back(dataset[indices[i]]);
        test_labels.push_back(labels[indices[i]]);
    }
    //残りをトレーニング用に格納
    for(int i=testsize;i<(int)dataset.size();i++){
        training_dataset.push_back(dataset[i]);
        training_labels.push_back(labels[indices[i]]);
    }
}

void Evaluation(vector<vector<double>> test_dataset, vector<int> test_labels){
    int TP=0,FP=0,FN=0,TN=0;

    for(int i=0;i<(int)test_dataset.size();i++){
        int predict=1;
        int actual=test_labels[i];
        if(predict==1&&actual==1){
            TP++;
        }else if(predict==1&&actual==0){
            FP++;
        }else if(predict==0&&actual==1){
            FN++;
        }else{
            TN++;
        }
    }
    int total = TP + FP + FN + TN;
    double accuracy=(double)(TP+TN)/total;
    double precision=(TP+FP>0)?(double)TP/(TP+FP):0.0;
    double recall=(TP + FN > 0)? (double)TP/(TP+FN):0.0;
    double f_score=(precision + recall > 0)? 2.0*precision * recall/(precision+recall) :0.0;

    cout << "Accuracy: "  << accuracy  << "\n";
    cout << "Precision: " << precision << "\n";
    cout << "Recall: "    << recall    << "\n";
    cout << "F-score: "   << f_score   << "\n";
    cout << "Confusion Matrix\n";
    cout << "TP: " << TP << "  FP: " << FP << "\n";
    cout << "FN: " << FN << "  TN: " << TN << "\n";
    

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
    Evaluation(test_dataset, test_labels);

}