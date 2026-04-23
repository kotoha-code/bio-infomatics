// ==============================================================
//  タンパク質可溶性予測 – 決定木 実装1〜4
//  コンパイル例: g++ -O2 -std=c++17 -o sol solubility_decision_tree.cpp
// ==============================================================

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

// ---- マクロ定義 ----
#define NUM_FEATURES 53    // 特徴量の次元数
#define NUM_SEQS     10000 // データセット総数

// ============================================================
//  TreeNode 構造体
//  (実装3 で定義・使用、実装4 でも共通利用)
// ============================================================
struct TreeNode {
    int    feature_id;    // 分岐に使う特徴量のインデックス
    double threshold;     // 分岐閾値（以下なら左、超えたら右）
    int    left_class_id; // 左に落ちたデータの予測ラベル
    int    right_class_id;// 右に落ちたデータの予測ラベル
};

// ============================================================
//  実装1(前半): LoadSolubilityFile
//    ファイルを読み込み、dataset と labels を埋める
//  ファイル形式:
//    1行目: タンパク質ID, 特徴量名×53, label
//    2行目〜: ID, 特徴量値×53, 0 or 1
// ============================================================
void LoadSolubilityFile(const string& filename,
                        vector<string>& feature_name,
                        vector<vector<double>>& dataset,
                        vector<int>& labels)
{
    ifstream ifs(filename);
    if (!ifs) {
        cerr << "[Error] Cannot open: " << filename << endl;
        return;
    }

    string line;

    // --- 1行目: ヘッダー ---
    getline(ifs, line);
    istringstream hss(line);
    string tok;
    vector<string> headers;
    while (getline(hss, tok, ','))
        headers.push_back(tok);
    // headers[0]=ID, headers[1..53]=特徴量名, headers[54]=label
    for (int i = 1; i <= NUM_FEATURES; i++)
        feature_name[i - 1] = headers[i];

    // --- 2行目以降: データ ---
    int row = 0;
    while (getline(ifs, line) && row < NUM_SEQS) {
        istringstream ss(line);
        getline(ss, tok, ','); // タンパク質IDを読み飛ばす

        for (int j = 0; j < NUM_FEATURES; j++) {
            getline(ss, tok, ',');
            dataset[row][j] = stod(tok);
        }
        getline(ss, tok, ','); // ラベル列
        labels[row] = stoi(tok);
        row++;
    }
}

// ============================================================
//  実装1(後半): DivideDataset
//    インデックス配列を作りシャッフル後、
//    先頭 test_ratio 割をテスト、残りをトレーニングに分割
// ============================================================
void DivideDataset(const vector<vector<double>>& dataset,
                   const vector<int>& labels,
                   vector<vector<double>>& training_dataset,
                   vector<int>& training_labels,
                   vector<vector<double>>& test_dataset,
                   vector<int>& test_labels,
                   double test_ratio)
{
    int n = (int)dataset.size();
    vector<int> indices(n);
    iota(indices.begin(), indices.end(), 0); // 0,1,2,...,n-1

    // seed を固定してシャッフル（実行ごとに結果が変わらないように）
    mt19937 rng(42);
    shuffle(indices.begin(), indices.end(), rng);

    int test_size = (int)(n * test_ratio); // 先頭 20% をテスト用

    for (int i = 0; i < test_size; i++) {
        test_dataset.push_back(dataset[indices[i]]);
        test_labels.push_back(labels[indices[i]]);
    }
    for (int i = test_size; i < n; i++) {
        training_dataset.push_back(dataset[indices[i]]);
        training_labels.push_back(labels[indices[i]]);
    }
}

// ============================================================
//  共通ユーティリティ: 混同行列の表示
// ============================================================
void PrintMetrics(int TP, int FP, int FN, int TN) {
    int total = TP + FP + FN + TN;
    double accuracy  = (double)(TP + TN) / total;
    double precision = (TP + FP > 0) ? (double)TP / (TP + FP) : 0.0;
    double recall    = (TP + FN > 0) ? (double)TP / (TP + FN) : 0.0;
    double f_score   = (precision + recall > 0)
                       ? 2.0 * precision * recall / (precision + recall)
                       : 0.0;

    cout << "Accuracy: "  << accuracy  << "\n";
    cout << "Precision: " << precision << "\n";
    cout << "Recall: "    << recall    << "\n";
    cout << "F-score: "   << f_score   << "\n";
    cout << "Confusion Matrix\n";
    cout << "TP: " << TP << "  FP: " << FP << "\n";
    cout << "FN: " << FN << "  TN: " << TN << "\n";
}

// ============================================================
//  実装2: Evaluation（全データをラベル1と予測する場合）
//    Recall=1 になるが Precision は低い
// ============================================================
void Evaluation(const vector<vector<double>>& test_dataset,
                const vector<int>& test_labels)
{
    int TP = 0, FP = 0, FN = 0, TN = 0;
    for (int i = 0; i < (int)test_dataset.size(); i++) {
        int predicted = 1; // 全て「可溶」と予測
        int actual    = test_labels[i];
        if      (predicted == 1 && actual == 1) TP++;
        else if (predicted == 1 && actual == 0) FP++;
        else if (predicted == 0 && actual == 1) FN++;
        else                                    TN++;
    }
    PrintMetrics(TP, FP, FN, TN);
}
// 期待出力: Accuracy=0.4075, Precision=0.4075, Recall=1,
//           F=0.579041, TP=815, FP=1185, FN=0, TN=0

// ============================================================
//  実装3 補助: ジニ不純度の計算
//    G = 2 * p * (1-p)  (p = ラベル1の割合)
//    最大値 0.5（p=0.5 のとき）、最小値 0（どちらかのみ）
// ============================================================
static double GiniImpurity(const vector<int>& lbls) {
    if (lbls.empty()) return 0.0;
    int cnt1 = 0;
    for (int l : lbls) cnt1 += l;
    double p = (double)cnt1 / (int)lbls.size();
    return 2.0 * p * (1.0 - p);
}

// 重み付きジニ不純度: (Nl/N)*Gl + (Nr/N)*Gr
static double WeightedGini(const vector<int>& L,
                            const vector<int>& R)
{
    int nl = (int)L.size(), nr = (int)R.size(), n = nl + nr;
    if (n == 0) return 0.0;
    return (double)nl / n * GiniImpurity(L) +
           (double)nr / n * GiniImpurity(R);
}

// ============================================================
//  実装3: TrainDecisionNode
//    全特徴量 × 1〜99パーセンタイル を閾値候補として探索し、
//    重み付きジニ不純度が最小になる (特徴量, 閾値) を選ぶ
// ============================================================
void TrainDecisionNode(const vector<vector<double>>& dataset,
                       const vector<int>& labels,
                       TreeNode& node)
{
    int n = (int)dataset.size();
    double best_gini = 1e18;

    for (int f = 0; f < NUM_FEATURES; f++) {
        // 特徴量 f の値を昇順ソート → パーセンタイル計算に使う
        vector<double> vals(n);
        for (int i = 0; i < n; i++) vals[i] = dataset[i][f];
        sort(vals.begin(), vals.end());

        // 1〜99 パーセンタイルを候補閾値として評価（計 99 個）
        for (int pct = 1; pct <= 99; pct++) {
            int idx = (int)((double)pct / 100.0 * (n - 1));
            double thr = vals[idx];

            vector<int> left_lbl, right_lbl;
            for (int i = 0; i < n; i++) {
                if (dataset[i][f] <= thr) left_lbl.push_back(labels[i]);
                else                      right_lbl.push_back(labels[i]);
            }

            double g = WeightedGini(left_lbl, right_lbl);
            if (g < best_gini) {
                best_gini = g;
                node.feature_id = f;
                node.threshold  = thr;

                // 左右それぞれで多数決ラベルを決める
                int l1 = 0; for (int l : left_lbl)  l1 += l;
                int r1 = 0; for (int l : right_lbl)  r1 += l;
                node.left_class_id  = (l1 * 2 >= (int)left_lbl.size())  ? 1 : 0;
                node.right_class_id = (r1 * 2 >= (int)right_lbl.size()) ? 1 : 0;
            }
        }
    }
}

// ============================================================
//  実装3: Evaluation（深さ1 TreeNode を使った予測）
// ============================================================
void Evaluation(const TreeNode& node,
                const vector<vector<double>>& test_dataset,
                const vector<int>& test_labels)
{
    int TP = 0, FP = 0, FN = 0, TN = 0;
    for (int i = 0; i < (int)test_dataset.size(); i++) {
        int predicted = (test_dataset[i][node.feature_id] <= node.threshold)
                        ? node.left_class_id
                        : node.right_class_id;
        int actual = test_labels[i];
        if      (predicted == 1 && actual == 1) TP++;
        else if (predicted == 1 && actual == 0) FP++;
        else if (predicted == 0 && actual == 1) FN++;
        else                                    TN++;
    }
    PrintMetrics(TP, FP, FN, TN);
}

// ============================================================
//  実装4: TrainDecisionTree（深さ2）
//    ノード0（親）→ TrainDecisionNode で学習
//    親の分割結果で左右データを分け、
//    ノード1（左の子）・ノード2（右の子）を個別に学習
// ============================================================
void TrainDecisionTree(const vector<vector<double>>& dataset,
                       const vector<int>& labels,
                       vector<TreeNode>& decision_tree)
{
    // --- 親ノード (decision_tree[0]) の学習 ---
    TrainDecisionNode(dataset, labels, decision_tree[0]);

    // --- 親ノードでデータを左右に分割 ---
    vector<vector<double>> left_ds, right_ds;
    vector<int>                 left_lbl, right_lbl;

    for (int i = 0; i < (int)dataset.size(); i++) {
        if (dataset[i][decision_tree[0].feature_id] <= decision_tree[0].threshold) {
            left_ds.push_back(dataset[i]);
            left_lbl.push_back(labels[i]);
        } else {
            right_ds.push_back(dataset[i]);
            right_lbl.push_back(labels[i]);
        }
    }

    // --- 左の子ノード (decision_tree[1]) の学習 ---
    TrainDecisionNode(left_ds,  left_lbl,  decision_tree[1]);
    // --- 右の子ノード (decision_tree[2]) の学習 ---
    TrainDecisionNode(right_ds, right_lbl, decision_tree[2]);
}

// ============================================================
//  実装4: Evaluation（深さ2 vector<TreeNode> を使った予測）
// ============================================================
void Evaluation(const vector<TreeNode>& decision_tree,
                const vector<vector<double>>& test_dataset,
                const vector<int>& test_labels)
{
    int TP = 0, FP = 0, FN = 0, TN = 0;
    for (int i = 0; i < (int)test_dataset.size(); i++) {
        // ルートノード[0] で左右を判定
        int predicted;
        if (test_dataset[i][decision_tree[0].feature_id] <= decision_tree[0].threshold) {
            // 左の子ノード[1] で最終予測
            const TreeNode& child = decision_tree[1];
            predicted = (test_dataset[i][child.feature_id] <= child.threshold)
                        ? child.left_class_id : child.right_class_id;
        } else {
            // 右の子ノード[2] で最終予測
            const TreeNode& child = decision_tree[2];
            predicted = (test_dataset[i][child.feature_id] <= child.threshold)
                        ? child.left_class_id : child.right_class_id;
        }

        int actual = test_labels[i];
        if      (predicted == 1 && actual == 1) TP++;
        else if (predicted == 1 && actual == 0) FP++;
        else if (predicted == 0 && actual == 1) FN++;
        else                                    TN++;
    }
    PrintMetrics(TP, FP, FN, TN);
}

// ============================================================
//  main
// ============================================================
int main(void) {

    // ---- 実装1: データ読み込みと分割 ----
    vector<string>          feature_name(NUM_FEATURES, "");
    vector<vector<double>>  dataset(NUM_SEQS, vector<double>(NUM_FEATURES, 0.0));
    vector<int>                  labels(NUM_SEQS);

    LoadSolubilityFile("protein_solubility_dataset.txt",
                       feature_name, dataset, labels);

    vector<vector<double>> training_dataset, test_dataset;
    vector<int>                 training_labels,  test_labels;
    double test_ratio = 0.2;

    DivideDataset(dataset, labels,
                  training_dataset, training_labels,
                  test_dataset,     test_labels,
                  test_ratio);

    // ---- 実装2: 全部ラベル1と予測した場合の評価 ----
    cout << "=== 実装2: 全部1予測 ===" << endl;
    Evaluation(test_dataset, test_labels);

    // ---- 実装3: 深さ1の決定木 ----
    cout << "\n=== 実装3: 深さ1の決定木 ===" << endl;
    TreeNode decision_tree;
    TrainDecisionNode(training_dataset, training_labels, decision_tree);
    Evaluation(decision_tree, test_dataset, test_labels);

    cout << "  使用特徴量ID: " << decision_tree.feature_id
              << "  閾値: "        << decision_tree.threshold
              << "  左ラベル: "    << decision_tree.left_class_id
              << "  右ラベル: "    << decision_tree.right_class_id << endl;

    // ---- 実装4: 深さ2の決定木 ----
    cout << "\n=== 実装4: 深さ2の決定木 ===" << endl;
    vector<TreeNode> decision_tree2(3);
    TrainDecisionTree(training_dataset, training_labels, decision_tree2);
    Evaluation(decision_tree2, test_dataset, test_labels);

    // 決定木の構造を表示（レポート用）
    cout << "\n--- 決定木の構造 ---\n";
    for (int k = 0; k < 3; k++) {
        cout << "ノード" << k
                  << " | feature=" << decision_tree2[k].feature_id
                  << " (" << feature_name[decision_tree2[k].feature_id] << ")"
                  << " | threshold=" << decision_tree2[k].threshold
                  << " | left="      << decision_tree2[k].left_class_id
                  << " | right="     << decision_tree2[k].right_class_id
                  << endl;
    }

    return 0;
}