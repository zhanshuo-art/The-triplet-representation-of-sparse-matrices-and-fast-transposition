#include <bits/stdc++.h>
using namespace std;

// ========== 配置 ==========
const string PATH_DeepDiveAI   = "E:\\InnovationDataset\\DeepInnovationAI\\DeepDiveAI.csv";
const string PATH_DeepPatentAl = "E:\\InnovationDataset\\DeepInnovationAI\\DeepPatentAI.csv";

const int DEEPDIVE_FIRST_N = 1000;

// 三元组结构体定义
struct Triple {
    int row;    // 行号
    int col;    // 列号
    int value;  // 非零元素值

    Triple() : row(0), col(0), value(0) {}
    Triple(int r, int c, int v) : row(r), col(c), value(v) {}
};

// 三元组顺序表结构
struct TSMatrix {
    vector<Triple> data;  // 三元组数组
    int rows;             // 矩阵行数
    int cols;             // 矩阵列数
    int nums;             // 非零元素个数

    TSMatrix() : rows(0), cols(0), nums(0) {}
    TSMatrix(int r, int c) : rows(r), cols(c), nums(0) {}
};

// 快速转置算法
TSMatrix fastTranspose(const TSMatrix& M) {
    TSMatrix T(M.cols, M.rows);  // 转置矩阵的行列互换
    T.nums = M.nums;

    if (M.nums == 0) return T;

    // num[col]: 统计矩阵M的第col列中非零元素的个数
    vector<int> num(M.cols, 0);
    for (int t = 0; t < M.nums; ++t) {
        num[M.data[t].col]++;
    }

    // cpot[col]: 指示M中第col列的第一个非零元素的位置
    vector<int> cpot(M.cols, 0);
    cpot[0] = 0;
    for (int col = 1; col < M.cols; ++col) {
        cpot[col] = cpot[col - 1] + num[col - 1];
    }

    // 按M.data的顺序依次将每个三元组转置后放入T.data的适当位置
    T.data.resize(M.nums);
    for (int p = 0; p < M.nums; ++p) {
        int col = M.data[p].col;
        int q = cpot[col];
        T.data[q].row = M.data[p].col;
        T.data[q].col = M.data[p].row;
        T.data[q].value = M.data[p].value;
        cpot[col]++;
    }

    return T;
}

// 输出三元组顺序表到CSV文件
void outputTSMatrix(const TSMatrix& M, const string& filename,
                    const vector<string>& rowNames,
                    const vector<string>& colNames) {
    ofstream fout(filename);
    fout << "row_index,row_name,col_index,col_name,value\n";

    for (int i = 0; i < M.nums; ++i) {
        fout << M.data[i].row << ","
             << "\"" << rowNames[M.data[i].row] << "\","
             << M.data[i].col << ","
             << colNames[M.data[i].col] << ","
             << M.data[i].value << "\n";
    }

    fout.close();
    cout << "三元组顺序表已保存: " << filename << " (共 " << M.nums << " 个非零元素)\n";
}

// 简单 CSV 解析：按逗号分割，去除首尾引号和空格
vector<string> parseCSV(const string& line) {
    vector<string> fields;
    stringstream ss(line);
    string field;
    while (getline(ss, field, ',')) {
        // 去除首尾空格
        size_t start = 0, end = field.size();
        while (start < end && isspace(field[start])) start++;
        while (end > start && isspace(field[end-1])) end--;
        field = field.substr(start, end - start);
        // 去除首尾引号
        if (field.size() >= 2 && field[0] == '"' && field.back() == '"') {
            field = field.substr(1, field.size() - 2);
        }
        fields.push_back(field);
    }
    return fields;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout.setf(ios::unitbuf); // 禁用输出缓冲，立即输出

    cout << "========== 开始构建 PaperTechMatrix ==========\n\n";

    // 1) 读取 DeepDiveAI.csv 前1000条
    cout << "读取 DeepDiveAI.csv...\n";
    ifstream finDeepDive(PATH_DeepDiveAI);
    if (!finDeepDive.is_open()) {
        cerr << "[ERROR] 无法打开: " << PATH_DeepDiveAI << endl;
        return 1;
    }

    string line;
    getline(finDeepDive, line); // 跳过表头

    vector<string> paperList;
    unordered_set<string> paperSet;
    int readCount = 0;

    while (readCount < DEEPDIVE_FIRST_N && getline(finDeepDive, line)) {
        vector<string> cols = parseCSV(line);
        if (cols.size() >= 2) {
            string paperID = cols[1];
            if (!paperID.empty() && paperSet.find(paperID) == paperSet.end()) {
                paperList.push_back(paperID);
                paperSet.insert(paperID);
            }
        }
        readCount++;
        if (readCount % 100 == 0) {
            cout << "已读取 " << readCount << " 行\n";
            cout.flush();
        }
    }
    finDeepDive.close();
    cout << "DeepDiveAI 读取完毕: " << readCount << " 行, 唯一论文数: " << paperList.size() << "\n\n";
    cout.flush();

    // 建立论文索引
    unordered_map<string,int> paperIndex;
    for (size_t i = 0; i < paperList.size(); ++i) {
        paperIndex[paperList[i]] = (int)i;
    }

    // 2) 读取 DeepPatentAl.csv
    cout << "读取 DeepPatentAI.csv...\n";

    unordered_map<string, vector<string>> paperToIPC;

    ifstream finPatent(PATH_DeepPatentAl);
    if (!finPatent.is_open()) {
        cerr << "[ERROR] 无法打开: " << PATH_DeepPatentAl << endl;
        return 1;
    }

    getline(finPatent, line); // 跳过表头

    int lineNo = 0;
    int sampleCount = 0;
    while (getline(finPatent, line)) {
        lineNo++;
        if (line.empty()) continue;

        vector<string> cols = parseCSV(line);
        if (cols.size() >= 4) {
            string paperID = cols[1];
            string ipc = cols[3];
            if (!ipc.empty() && !paperID.empty()) {
                size_t pos = ipc.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
                if (pos != string::npos && pos + 4 <= ipc.size()) {
                    string ipc4 = ipc.substr(pos, 4);
                    paperToIPC[paperID].push_back(ipc4);

                    // 输出前3个样例
                    if (sampleCount < 3) {
                        cout << "样例" << (sampleCount+1) << " - 原始IPC: [" << ipc << "] -> 解析结果: [" << ipc4 << "]\n";
                        sampleCount++;
                    }
                }
            }
        }
        if (lineNo % 1000 == 0) {
            cout << "已读取 " << lineNo << " 行\n";
            cout.flush();
        }
    }
    finPatent.close();
    cout << "DeepPatentAI 读取完毕: " << lineNo << " 行\n\n";
    cout.flush();

    // 3) 构建矩阵
    cout << "构建论文-技术矩阵...\n";
    unordered_map<string,int> ipcIndex;
    vector< unordered_map<int,int> > matrix(paperList.size());

    auto getIpcIndex = [&](const string &ipc4) -> int {
        auto it = ipcIndex.find(ipc4);
        if (it != ipcIndex.end()) return it->second;
        int idx = (int)ipcIndex.size();
        ipcIndex[ipc4] = idx;
        return idx;
    };

    for (size_t ii = 0; ii < paperList.size(); ++ii) {
        const string &paperID = paperList[ii];
        auto it = paperToIPC.find(paperID);
        if (it != paperToIPC.end()) {
            for (const string &ipc4 : it->second) {
                int j = getIpcIndex(ipc4);
                matrix[ii][j]++;
            }
        }
        if (ii % 100 == 0) {
            cout << "已处理 " << ii << "/" << paperList.size() << " 个论文\n";
            cout.flush();
        }
    }

    vector<string> ipcList(ipcIndex.size());
    for (auto &kv : ipcIndex) {
        ipcList[kv.second] = kv.first;
    }
    cout << "矩阵构建完成, IPC分类数: " << ipcList.size() << "\n\n";
    cout.flush();

    // 4) 输出密集矩阵
    cout << "输出密集矩阵...\n";
    cout.flush();
    string outDense = "C:\\C++ code\\InstTechMatrix.csv";
    ofstream foutDense(outDense);

    foutDense << "paper_id";
    for (auto &ipc : ipcList) foutDense << "," << ipc;
    foutDense << "\n";
    foutDense.flush();

    for (size_t ii = 0; ii < paperList.size(); ++ii) {
        foutDense << "\"" << paperList[ii] << "\"";
        for (size_t j = 0; j < ipcList.size(); ++j) {
            int val = 0;
            auto it = matrix[ii].find((int)j);
            if (it != matrix[ii].end()) val = it->second;
            foutDense << "," << val;
        }
        foutDense << "\n";
        if ((ii + 1) % 50 == 0) {
            cout << "已输出 " << (ii + 1) << "/" << paperList.size() << " 个论文\n";
            cout.flush();
            foutDense.flush();
        }
    }
    foutDense.close();
    cout << "密集矩阵已保存: " << outDense << "\n\n";
    cout.flush();

    // 5) 将矩阵转换为三元组顺序表
    cout << "\n========== 开始三元组转换与转置 ==========\n\n";
    cout << "将矩阵转换为三元组顺序表...\n";

    TSMatrix A(paperList.size(), ipcList.size());

    // 遍历矩阵，提取所有非零元素
    for (size_t i = 0; i < paperList.size(); ++i) {
        for (auto &p : matrix[i]) {
            int j = p.first;
            int val = p.second;
            if (val != 0) {
                A.data.push_back(Triple(i, j, val));
                A.nums++;
            }
        }
        if ((i + 1) % 100 == 0) {
            cout << "已转换 " << (i + 1) << "/" << paperList.size() << " 个论文\n";
            cout.flush();
        }
    }

    cout << "三元组顺序表构建完成\n";
    cout << "矩阵规模: " << A.rows << " × " << A.cols << "\n";
    cout << "非零元素个数: " << A.nums << "\n\n";
    cout.flush();

    // 输出原始三元组顺序表
    cout << "输出原始矩阵的三元组顺序表...\n";
    outputTSMatrix(A, "C:\\C++ code\\InstTechMatrix_triple.csv", paperList, ipcList);
    cout.flush();

    // 6) 快速转置算法
    cout << "\n执行快速转置算法...\n";

    TSMatrix AT = fastTranspose(A);

    cout << "转置矩阵规模: " << AT.rows << " × " << AT.cols << "\n";
    cout << "非零元素个数: " << AT.nums << "\n\n";
    cout.flush();

    // 输出转置后的三元组顺序表
    cout << "输出转置矩阵的三元组顺序表...\n";
    outputTSMatrix(AT, "C:\\C++ code\\InstTechMatrix_transposed.csv", ipcList, paperList);
    cout.flush();

    // ================================================================
    // 7) 计算机构合作潜力矩阵 (A × A^T)
    // ================================================================
    cout << "\n========== 计算机构合作潜力矩阵 ==========\n\n";
    cout << "计算矩阵乘法 A × A^T...\n";
    cout << "结果矩阵规模: " << paperList.size() << " × " << paperList.size() << "\n";
    cout.flush();

    auto mult_start = chrono::high_resolution_clock::now();

    // 使用稀疏矩阵乘法计算 A × A^T
    // 结果矩阵: collaboration[i][j] 表示论文i和论文j的技术重叠度
    vector<vector<int>> collaboration(paperList.size(), vector<int>(paperList.size(), 0));

    // 对于每对论文 (i, j)，计算它们在所有技术分类上的乘积和
    for (size_t i = 0; i < paperList.size(); ++i) {
        for (size_t j = 0; j < paperList.size(); ++j) {
            int sum = 0;
            // 遍历所有技术分类
            for (size_t k = 0; k < ipcList.size(); ++k) {
                int val_i = 0, val_j = 0;

                // 获取 A[i][k]
                auto it_i = matrix[i].find((int)k);
                if (it_i != matrix[i].end()) val_i = it_i->second;

                // 获取 A[j][k] (即 A^T[k][j])
                auto it_j = matrix[j].find((int)k);
                if (it_j != matrix[j].end()) val_j = it_j->second;

                sum += val_i * val_j;
            }
            collaboration[i][j] = sum;
        }

        if ((i + 1) % 50 == 0) {
            cout << "已计算 " << (i + 1) << "/" << paperList.size() << " 行\n";
            cout.flush();
        }
    }

    auto mult_end = chrono::high_resolution_clock::now();
    auto mult_duration = chrono::duration_cast<chrono::milliseconds>(mult_end - mult_start);

    cout << "矩阵乘法完成，耗时: " << mult_duration.count() << " 毫秒\n\n";
    cout.flush();

    // 输出机构合作潜力矩阵
    cout << "输出机构合作潜力矩阵...\n";
    string outCollaboration = "C:\\C++ code\\InstCollaborationMatrix.csv";
    ofstream foutCollab(outCollaboration);

    // 输出表头（第一行是所有论文ID）
    foutCollab << "paper_id";
    for (size_t j = 0; j < paperList.size(); ++j) {
        foutCollab << ",\"" << paperList[j] << "\"";
    }
    foutCollab << "\n";
    foutCollab.flush();

    // 输出每一行
    for (size_t i = 0; i < paperList.size(); ++i) {
        foutCollab << "\"" << paperList[i] << "\"";
        for (size_t j = 0; j < paperList.size(); ++j) {
            foutCollab << "," << collaboration[i][j];
        }
        foutCollab << "\n";

        if ((i + 1) % 50 == 0) {
            cout << "已输出 " << (i + 1) << "/" << paperList.size() << " 行\n";
            cout.flush();
            foutCollab.flush();
        }
    }

    foutCollab.close();
    cout << "机构合作潜力矩阵已保存: " << outCollaboration << "\n\n";
    cout.flush();

    cout << "\n========== 所有任务完成 ==========\n";
    cout << "生成的文件:\n";
    cout << "1. InstTechMatrix.csv - 密集矩阵表示 (论文×技术)\n";
    cout << "2. InstTechMatrix_triple.csv - 原始矩阵的三元组顺序表\n";
    cout << "3. InstTechMatrix_transposed.csv - 转置矩阵的三元组顺序表\n";
    cout << "4. InstCollaborationMatrix.csv - 机构合作潜力矩阵 (论文×论文)\n\n";

    return 0;
}

