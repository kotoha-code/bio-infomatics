// ステップ1: マルチプルアラインメントから対数オッズスコア行列を計算する


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>

// 塩基 → インデックス  A=0, C=1, G=2, T=3
int baseToIndex(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:             return -1;
    }
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "使い方: " << argv[0] << " <motifファイル>" << std::endl;
        return 1;
    }
    std::string filename = argv[1];

    // ----------------------------------------------------------
    // ファイル読み込み: 1行1配列
    // ----------------------------------------------------------
    std::ifstream fin(filename);
    if (!fin) {
        std::cerr << "エラー: " << filename << " を開けません" << std::endl;
        return 1;
    }

    std::vector<std::string> seqs;
    std::string line;
    while (std::getline(fin, line)) {
        while (!line.empty() && (line.back()=='\r' || line.back()==' '))
            line.pop_back();
        if (!line.empty()) seqs.push_back(line);
    }

    if (seqs.empty()) {
        std::cerr << "エラー: 配列が読み込めませんでした" << std::endl;
        return 1;
    }

    int N = (int)seqs.size();   // 配列本数
    int L = (int)seqs[0].size();// モチーフ長

    std::cout << "=== 入力データ ===" << std::endl;
    std::cout << "ファイル   : " << filename << std::endl;
    std::cout << "配列本数   : " << N << std::endl;
    std::cout << "モチーフ長 : " << L << std::endl;
    std::cout << std::endl;

    // ----------------------------------------------------------
    // 2次元配列を vector<vector<int/double>> で確保
    // freq[位置][塩基],  prob[位置][塩基],  score[位置][塩基]
    // ----------------------------------------------------------
    // ゼロ初期化
    std::vector<std::vector<int>>    freq (L, std::vector<int>   (4, 0));
    std::vector<std::vector<double>> prob (L, std::vector<double>(4, 0.0));
    std::vector<std::vector<double>> score(L, std::vector<double>(4, 0.0));

    // ----------------------------------------------------------
    // ステップ1: 各列の塩基出現回数を数える（頻度表）
    // ----------------------------------------------------------
    for (const std::string& seq : seqs) {
        for (int j = 0; j < (int)seq.size() && j < L; j++) {
            int b = baseToIndex(seq[j]);
            if (b >= 0) freq[j][b]++;
        }
    }

    const char* BASE_NAMES = "ACGT";

    std::cout << "=== 頻度表（+1前） ===" << std::endl;
    std::cout << "   ";
    for (int j = 0; j < L; j++)
        std::cout << std::setw(5) << "(" << j+1 << ")";
    std::cout << std::endl;
    for (int b = 0; b < 4; b++) {
        std::cout << BASE_NAMES[b] << "  ";
        for (int j = 0; j < L; j++)
            std::cout << std::setw(6) << freq[j][b];
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // ----------------------------------------------------------
    // ステップ2: 擬似頻度として各要素に +1 する
    //   → freq[j][b] が 0 でも log(0) にならないようにする
    // ----------------------------------------------------------
    for (int j = 0; j < L; j++)
        for (int b = 0; b < 4; b++)
            freq[j][b] += 1;

    std::cout << "=== 頻度表（+1後） ===" << std::endl;
    std::cout << "   ";
    for (int j = 0; j < L; j++)
        std::cout << std::setw(5) << "(" << j+1 << ")";
    std::cout << std::endl;
    for (int b = 0; b < 4; b++) {
        std::cout << BASE_NAMES[b] << "  ";
        for (int j = 0; j < L; j++)
            std::cout << std::setw(6) << freq[j][b];
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // ----------------------------------------------------------
    // ステップ3: 出現確率 p_i(x) を計算する
    //   p_i(x) = freq[i][x] / 列の合計
    // ----------------------------------------------------------
    for (int j = 0; j < L; j++) {
        int col_total = freq[j][0] + freq[j][1] + freq[j][2] + freq[j][3];
        for (int b = 0; b < 4; b++)
            prob[j][b] = (double)freq[j][b] / col_total;
    }

    // ----------------------------------------------------------
    // ステップ4: バックグラウンド確率 q(x) を計算する
    //   全結合配列中の塩基出現頻度から算出
    // ----------------------------------------------------------
    long long bg_count[4] = {0, 0, 0, 0};
    long long bg_total = 0;
    for (const std::string& seq : seqs) {
        for (char c : seq) {
            int b = baseToIndex(c);
            if (b >= 0) { bg_count[b]++; bg_total++; }
        }
    }

    double q[4];
    for (int b = 0; b < 4; b++)
        q[b] = (double)bg_count[b] / bg_total;

    std::cout << "=== バックグラウンド確率 q(x) ===" << std::endl;
    for (int b = 0; b < 4; b++)
        std::cout << "  q(" << BASE_NAMES[b] << ") = "
                  << std::fixed << std::setprecision(4) << q[b] << std::endl;
    std::cout << std::endl;

    // ----------------------------------------------------------
    // ステップ5: 対数オッズスコア s_i(x) = log2( p_i(x) / q(x) )
    // ----------------------------------------------------------
    for (int j = 0; j < L; j++)
        for (int b = 0; b < 4; b++)
            score[j][b] = log2(prob[j][b] / q[b]);

    std::cout << "=== 対数オッズスコア行列  s_i(x) = log2( p_i(x) / q(x) ) ===" << std::endl;
    std::cout << "   ";
    for (int j = 0; j < L; j++)
        std::cout << std::setw(6) << "(" << j+1 << ")";
    std::cout << std::endl;
    for (int b = 0; b < 4; b++) {
        std::cout << BASE_NAMES[b] << "  ";
        for (int j = 0; j < L; j++)
            std::cout << std::fixed << std::setprecision(2)
                      << std::setw(7) << score[j][b];
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // ----------------------------------------------------------
    // スコア行列をファイルに保存（ステップ2で読み込む）
    // フォーマット:
    //   1行目: モチーフ長 L
    //   次のL行: 位置ごとに A C G T のスコアをスペース区切り
    //   最終行: A C G T のバックグラウンド確率
    // ----------------------------------------------------------
    std::string out_file = filename + ".score";
    std::ofstream fout(out_file);
    fout << L << std::endl;
    for (int j = 0; j < L; j++) {
        for (int b = 0; b < 4; b++) {
            if (b > 0) fout << " ";
            fout << std::fixed << std::setprecision(6) << score[j][b];
        }
        fout << std::endl;
    }
    // バックグラウンド確率を最終行に保存
    for (int b = 0; b < 4; b++) {
        if (b > 0) fout << " ";
        fout << std::fixed << std::setprecision(6) << q[b];
    }
    fout << std::endl;

    std::cout << "スコア行列を " << out_file << " に保存しました" << std::endl;
    std::cout << "（このファイルをステップ2で使います）" << std::endl;

    return 0;
}
