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
    //cout<<headers.size()<<endl;
    //headerの中でid名以外を格納
    for(int i=1; i<=NUM_FEATURES; i++){
        feature_name[i-1]=headers[i];
    }

    int row=0;
    while(getline(ifs,line)&&row<NUM_SEQS){
        istringstream ss(line); //読んだ行をstreamに変換
        getline(ss,temp,' ');  //protein idは使わないから捨てる

        for(int i=0;i<NUM_FEATURES;i++){
            getline(ss,temp,' ');
            dataset[row][i]=stod(temp); //stod=文字列から少数に変換する
        }
        getline(ss,temp,' ');
        labels[row]=stoi(temp); //stoi=文字列から整数に変換する
        row++;
    }
}

void DivideDataset(vector<vector<double>>& dataset,vector<int>& labels, vector<vector<double>>& training_dataset,vector<int>& training_labels,vector<vector<double>>& test_dataset,vector<int>& test_labels,double test_ratio){
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
        training_dataset.push_back(dataset[indices[i]]);
        training_labels.push_back(labels[indices[i]]);
    }
}

void Evaluation(TreeNode node,vector<vector<double>> test_dataset, vector<int> test_labels){
    int TP=0,FP=0,FN=0,TN=0;
    for(int i=0;i<(int)test_dataset.size();i++){
        //int predict=1;
        int predict;
        if(test_dataset[i][node.feature_id] <= node.threshold){
            predict=node.left_class_id;
        }else{
            predict=node.right_class_id;
        }
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

    int total= TP+FP+FN+TN;
    double accuracy=(double)(TP+TN)/total;
    double precision=(TP+FP>0)?(double)TP/(TP+FP):0.0;
    double recall=(TP + FN > 0)? (double)TP/(TP+FN):0.0;
    double f_score=(precision + recall > 0)? 2.0*precision*recall/(precision+recall) :0.0;

    cout<<"Accuracy: "<< accuracy<<endl;
    cout<<"Precision: "<< precision<<endl;
    cout<<"Recall: "<< recall<<endl;
    cout<<"F-score: "<< f_score<<endl;
    cout<<"Confusion Matrix\n";
    cout<<"TP: "<<TP<<"  FP: "<<FP<<"\n";
    cout<<"FN: "<<FN<<"  TN: "<<TN<<"\n";
    
}

//ジニ不純度
double GiniImpurity(const vector<int>& x){
    int count=0;
    for (int i=0; i<x.size(); i++){
        count += x[i];
    }
    double p = (double)count / (int)x.size();
    return 2.0*p*(1.0-p);
}

// 重み付きジニ不純度
double WeightedGini(const vector<int>& L,const vector<int>& R){
    int nl=(int)L.size();
    int nr=(int)R.size();
    int n=nl+nr;
    if (n==0){
        return 0.0;
    }
    return (double)nl/n * GiniImpurity(L)+(double)nr/n * GiniImpurity(R);
}

void TrainDecisionNode(const vector<vector<double>>& dataset,const vector<int>& labels,TreeNode& node){
    int n=dataset.size();
    double best_gini =1e18; //最小値を求めたいので初期値を大きく設定

    for (int f = 0; f<NUM_FEATURES; f++){
        // 特徴量 f の値を昇順ソート → パーセンタイル計算に使う
        vector<double> vals(n);
        for (int i=0; i<n; i++) vals[i] = dataset[i][f];
        sort(vals.begin(), vals.end());

        // 1〜99 パーセンタイルを候補閾値として評価（計 99 個）
        for (int pct=1; pct<100; pct++) {
            int idx = (int)(n*(double)pct/100.0);
            double thr = vals[idx];
            vector<int> left;
            vector<int> right;
            for (int i=0; i<n; i++){
                if (dataset[i][f] <= thr){
                    left.push_back(labels[i]);
                }else{
                    right.push_back(labels[i]);
                }                     
            }

            double g = WeightedGini(left, right);

            if (g < best_gini) {
                best_gini = g;
                node.feature_id = f;
                node.threshold  = thr;
                node.left_class_id=0;
                node.right_class_id=0;

                // 左右それぞれで多数決ラベルを決める
                int l1=count(left.begin(),left.end(),1);
                if(l1>=left.size()/2){
                    node.left_class_id=1;
                }
                int r1=count(right.begin(),right.end(),1);
                if(r1>=right.size()/2){
                    node.right_class_id=1;
                }
            }
        }
    }
}

void Evaluation(vector<TreeNode>& node,vector<vector<double>> test_dataset, vector<int> test_labels){
    int TP=0,FP=0,FN=0,TN=0;
    for(int i=0;i<(int)test_dataset.size();i++){
        //int predict=1;
        int predict;
        
        if(test_dataset[i][node[0].feature_id] <= node[0].threshold){

            if(test_dataset[i][node[1].feature_id] <= node[1].threshold){
                predict=node[1].left_class_id;
            }else{
                predict=node[1].right_class_id;
            }
        }else{

            if(test_dataset[i][node[2].feature_id] <= node[2].threshold){
                predict=node[2].left_class_id;
            }else{
                predict=node[2].right_class_id;
            }
        } 
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

    int total= TP+FP+FN+TN;
    double accuracy=(double)(TP+TN)/total;
    double precision=(TP+FP>0)?(double)TP/(TP+FP):0.0;
    double recall=(TP + FN > 0)? (double)TP/(TP+FN):0.0;
    double f_score=(precision + recall > 0)? 2.0*precision*recall/(precision+recall) :0.0;

    cout<<"Accuracy: "<< accuracy<<endl;
    cout<<"Precision: "<< precision<<endl;
    cout<<"Recall: "<< recall<<endl;
    cout<<"F-score: "<< f_score<<endl;
    cout<<"Confusion Matrix\n";
    cout<<"TP: "<<TP<<"  FP: "<<FP<<"\n";
    cout<<"FN: "<<FN<<"  TN: "<<TN<<"\n";
    
}

void TrainDecisionTree(const vector<vector<double>>& dataset,const vector<int>& labels,vector<TreeNode>& tree){
    TrainDecisionNode(dataset,labels,tree[0]);

    vector<vector<double>>left_dataset,right_dataset;
    vector<int>left_label,right_label;

    for(int i=0;i<dataset.size();i++){
        if(dataset[i][tree[0].feature_id] <= tree[0].threshold){
            left_dataset.push_back(dataset[i]);
            left_label.push_back(labels[i]);
        }else{
            right_dataset.push_back(dataset[i]);
            right_label.push_back(labels[i]);
        }
    }
    
    TrainDecisionNode(left_dataset,left_label,tree[1]);
    TrainDecisionNode(right_dataset,right_label,tree[2]);
}

int main(void){
    /*TreeNode decision_tree1;
    TreeNode decision_tree2;
    TreeNode decision_tree3;*/

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
   
    vector<TreeNode> decision_tree(3);

    TrainDecisionTree(training_dataset,training_labels,decision_tree);
    cout<<decision_tree.size()<<endl;

    Evaluation(decision_tree,test_dataset,test_labels);

    cout<<feature_name[decision_tree[0].feature_id]<<" : "<<decision_tree[0].feature_id<<" "<<decision_tree[0].threshold<<" "<<decision_tree[0].left_class_id<<" "<<decision_tree[0].right_class_id<<endl;
    cout<<feature_name[decision_tree[1].feature_id]<<" : "<< decision_tree[1].feature_id<<" "<<decision_tree[1].threshold<<" "<<decision_tree[1].left_class_id<<" "<<decision_tree[1].right_class_id<<endl;
    cout<<feature_name[decision_tree[2].feature_id]<<" : "<< decision_tree[2].feature_id<<" "<<decision_tree[2].threshold<<" "<<decision_tree[2].left_class_id<<" "<<decision_tree[2].right_class_id<<endl;

    return 0;

}