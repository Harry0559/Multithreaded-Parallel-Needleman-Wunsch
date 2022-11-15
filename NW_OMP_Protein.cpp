#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <omp.h>
#include <time.h>
#include<iomanip>
#include<string.h>

using namespace std;

void fillBlosum62(vector<vector<int>> &blosum62);

class myMatrix
{
private:
    string A;
    string B;
    int lenA;
    int lenB;
    int threadNums;
    string aminoAcid;
    vector<vector<int>> blosum62;
    vector<vector<int>> matrix;

public:
    myMatrix(string &x, string &y, int thread)
    {
        A = x;
        B = y;
        lenA = A.size();
        lenB = B.size();
        threadNums = thread;
        fillBlosum62(blosum62); //use BLOSUM62 as the substitution matrix
        aminoAcid = "*ARNDCQEGHILKMFPSTWYVBZX"; //record the order of amino acid in BLOSUM62
        for (int i = 0; i <= lenA; i++)
        {
            matrix.emplace_back(vector<int>(lenB + 1, 0));
        }
    }
    void initialize()
    {
        omp_set_num_threads(threadNums);
        #pragma omp parallel
        {
        #pragma omp for
            for (int i = 0; i <= lenA; i++)
            {
                matrix[i][0] = -4 * i;
            }
        #pragma omp for
            for (int j = 0; j <= lenB; j++)
            {
                matrix[0][j] = -4 * j;
            }
        }
    }
    int getIndex(char x) {
        int ans = 0;
        for (int i = 1; i < 24; i++) {
            if (aminoAcid[i] == x) {
                ans = i;
                break;
            }
        }
        return ans;
    }
    void fillElement(int i, int j)
    {
        if (i > lenA || j > lenB || i == 0 || j == 0) return; // comfirm that i and j are valid indices
        int S = blosum62[getIndex(A[i-1])][getIndex(B[j-1])];
        int match = matrix[i - 1][j - 1] + S;
        int del = matrix[i - 1][j] - 4;
        int ins = matrix[i][j - 1] - 4;
        matrix[i][j] = max(match, max(del, ins));
    }
    void fillMatrix()
    {
        int cycle = ceil((float)lenA / threadNums); // separation of the matrix
        for (int k = 0; k < cycle; k++)
        {
            int start = k * threadNums + 1;
            int thread = (k == cycle - 1) ? (lenA - start + 1) : threadNums;
            int end = start + thread - 1;
            // sequential part 1
            #pragma omp master
            {
                for (int i = 1; i < thread; i++)
                {
                    for (int j = 1; j < 1 + thread - i; j++)
                    {
                        fillElement(start + i - 1, j);
                    }
                }
            }
            // parallel part
            for (int count = 0; count < lenB - (thread - 1); count++) {
                omp_set_num_threads(thread);
                #pragma omp parallel for
                for (int i = start; i <= end; i++) {
                    fillElement(i, thread - (i - start) + count);
                }
            }
            // sequential part 2
            #pragma omp master
            {
                for (int i = 1; i < thread; i++)
                {
                    for (int j = lenB - (i - 1); j <= lenB; j++)
                    {
                        fillElement(start + i, j);
                    }
                }
            }
        }
    }
    void backtracking(string &alignmentA, string &alignmentB)
    {
        int i = lenA;
        int j = lenB;
        while (i > 0 && j > 0)
        {
            int score = matrix[i][j];
            int scoreDiag = matrix[i - 1][j - 1];
            int scoreUp = matrix[i][j - 1];
            int scoreLeft = matrix[i - 1][j];
            int S = blosum62[getIndex(A[i-1])][getIndex(B[j-1])];
            if (score == scoreDiag + S)
            {
                alignmentA = A[i - 1] + alignmentA;
                alignmentB = B[j - 1] + alignmentB;
                i--;
                j--;
            }
            else if (score == scoreLeft - 4)
            {
                alignmentA = A[i - 1] + alignmentA;
                alignmentB = "-" + alignmentB;
                i--;
            }
            else
            {
                alignmentA = "-" + alignmentA;
                alignmentB = B[j - 1] + alignmentB;
                j--;
            }
        }
        while (i > 0)
        {
            alignmentA = A[i - 1] + alignmentA;
            alignmentB = "-" + alignmentB;
            i--;
        }
        while (j > 0)
        {
            alignmentA = "-" + alignmentA;
            alignmentB = B[j - 1] + alignmentB;
            j--;
        }
    }
    void showMatrix()
    {
        string col = "--" + B;
        string row = "-" + A;
        for (int i = 0; i < col.size(); i++) cout << setw(5) <<col[i];
        cout << endl;
        for (int i = 0; i < row.size(); i++)
        {
            cout << setw(5) <<row[i];
            for (int j = 0; j < lenB + 1; j++) cout << setw(5) <<matrix[i][j];
            cout << endl;
        }
    }
};

vector<string> readSequences(string &pathA, string &pathB)
{
    vector<string> ans(2, string(0, ' '));
    string temp;
    ifstream fileA(pathA.c_str());
    getline(fileA, temp);
    while (!fileA.eof()) {
        getline(fileA, temp);
        ans[0].append(temp);
    }
    fileA.close();
    ifstream fileB(pathB.c_str());
    getline(fileB, temp);
    while (!fileB.eof()) {
        getline(fileB, temp);
        ans[1].append(temp);
    }
    fileB.close();
    return ans;
}


int main(int argc, char* argv[])
{
    int mode = 0;
    if (argc == 1 || (argc == 2 && strcmp(argv[1],"--interactive") == 0)) mode = 0;
    else if (argc == 2 && strcmp(argv[1],"--file") == 0) mode = 1;
    else if (argc == 2 && strcmp(argv[1],"--help") == 0) mode = 2;
    else if (argc == 2 && strcmp(argv[1],"--demo") == 0) mode = 3;
    else {
        cout << "The command is wrong. The correct usage is shown below:\n" << endl;
        mode = 2;
    }

    if (mode == 0) //interactive mode
    {
        vector<string> setA = {"data/3500.fasta","data/6995.fasta","data/13100.fasta","data/34350.fasta"};
        vector<string> setB = {"data/3503.fasta","data/7031.fasta","data/14507.fasta","data/35213.fasta"};
        int indexA;
        int indexB;
        cout << "length of sequence A(0:3500, 1:6995, 2:13100, 3:34350): ";
        cin >> indexA;
        cout << "length of sequence B(0:3503, 1:7031, 2:14507, 3:35213): ";
        cin >> indexB;
        cout << "\nThe sequence data is being read from the file..." << endl;
        vector<string> sequences = readSequences(setA[indexA], setB[indexB]);
        cout << "Reading has finished." << endl << endl;
        int threads;
        cout << "number of threads: ";
        cin >> threads;
        cout << "\nstart to time..." << endl;
        clock_t start = clock();
        myMatrix dp(sequences[0], sequences[1], threads);
        dp.initialize();
        dp.fillMatrix();
        string alignmentA;
        string alignmentB;
        dp.backtracking(alignmentA, alignmentB);
        clock_t end = clock();
        cout << "Execution time: " << (end - start) << "ms" << endl;
        system("pause");
    }
    else if (mode == 1) //file mode
    {
        vector<string> setA = {"data/3500.fasta","data/6995.fasta"};
        vector<string> setB = {"data/3503.fasta","data/7031.fasta"};
        vector<string> nameA = {"3500","6995"};
        vector<string> nameB = {"3503","7031"};
        vector<int> threads = {2,4,8,16};
        vector<string> sequences;
        clock_t start;
        clock_t end;
        clock_t execution_time;
        clock_t two_thread_time;
        string name;
        ofstream ofile;
        float speedUp;
        cout << setw(15) << "SeqA:SeqB" << setw(22) << "Number of threads" << setw(23) << "Execution time(ms)" << setw(12) << "Speedup" << endl;
        for (int i = 0; i < 2; i++) {
            sequences = readSequences(setA[i], setB[i]);
            for (int j = 0; j < 4; j++) {
                start = clock();
                myMatrix dp(sequences[0], sequences[1], threads[j]);
                dp.initialize();
                dp.fillMatrix();
                string alignmentA;
                string alignmentB;
                dp.backtracking(alignmentA, alignmentB);
                end = clock();
                execution_time = end - start;
                if (j == 0) two_thread_time = execution_time;
                speedUp = (float)two_thread_time / execution_time;
                name = "output/" + nameA[i] + "_" + nameB[i] + "_" + to_string(threads[j]) + ".txt";
                ofile.open(name.c_str());
                ofile << execution_time << endl;
                ofile << setprecision(3) << speedUp;
                ofile.close();
                cout << setw(15) << nameA[i]+":"+nameB[i] << setw(22) << threads[j] << setw(23) << execution_time << setw(12) << setprecision(3) << speedUp << endl;
            }
            if (i == 1) continue;
            sequences = readSequences(setA[i+1], setB[i]);
            for (int j = 0; j < 4; j++) {
                start = clock();
                myMatrix dp(sequences[0], sequences[1], threads[j]);
                dp.initialize();
                dp.fillMatrix();
                string alignmentA;
                string alignmentB;
                dp.backtracking(alignmentA, alignmentB);
                end = clock();
                execution_time = end - start;
                if (j == 0) two_thread_time = execution_time;
                speedUp = (float)two_thread_time / execution_time;
                name = "output/" + nameA[i+1] + "_" + nameB[i] + "_" + to_string(threads[j]) + ".txt";
                ofile.open(name.c_str());
                ofile << execution_time << endl;
                ofile << setprecision(3) << speedUp;
                ofile.close();
                cout << setw(15) << nameA[i+1]+":"+nameB[i] << setw(22) << threads[j] << setw(23) << execution_time << setw(12) << setprecision(3) << speedUp << endl;
            }
        }
        system("pause");
    }
    else if (mode == 2) //help mode
    {
        cout << "Usage:" << endl;
        cout << setw(40) << setiosflags(ios::left) << "./NW_OMP_Protein.exe" << "Interactive mode by default" << endl;
        cout << setw(40) << setiosflags(ios::left) << "./NW_OMP_Protein.exe --interactive" << "Interactive mode by explicit command" << endl;
        cout << setw(40) << setiosflags(ios::left) << "./NW_OMP_Protein.exe --file" << "File mode" << endl;
        cout << setw(40) << setiosflags(ios::left) << "./NW_OMP_Protein.exe --help" << "Show the correct usage of the command" << endl;
        cout << setw(40) << setiosflags(ios::left) << "./NW_OMP_Protein.exe --demo" << "Demonstrate the matrix to verify the correctness" << endl;
    }
    else if (mode == 3) //demonstration mode
    {
        string pathA = "data/seqA.fasta";
        string pathB = "data/seqB.fasta";
        vector<string> sequences = readSequences(pathA,pathB);
        int threads;
        cout << "Number of threads: ";
        cin >> threads;
        myMatrix dp(sequences[0], sequences[1], threads);
        dp.initialize();
        dp.fillMatrix();
        dp.showMatrix();
    }
}


void fillBlosum62(vector<vector<int>> &blosum62)
{
    for (int i = 0; i < 24; i++)
    {
        blosum62.emplace_back(vector<int>(24, -4));
    }
    blosum62[0][0] = 1;
    blosum62[1][1] = 4;
    blosum62[1][2] = -1;
    blosum62[1][3] = -2;
    blosum62[1][4] = -2;
    blosum62[1][5] = 0;
    blosum62[1][6] = -1;
    blosum62[1][7] = -1;
    blosum62[1][8] = 0;
    blosum62[1][9] = -2;
    blosum62[1][10] = -1;
    blosum62[1][11] = -1;
    blosum62[1][12] = -1;
    blosum62[1][13] = -1;
    blosum62[1][14] = -2;
    blosum62[1][15] = -1;
    blosum62[1][16] = 1;
    blosum62[1][17] = 0;
    blosum62[1][18] = -3;
    blosum62[1][19] = -2;
    blosum62[1][20] = 0;
    blosum62[1][21] = -2;
    blosum62[1][22] = -1;
    blosum62[1][23] = 0;
    blosum62[2][1] = -1;
    blosum62[2][2] = 5;
    blosum62[2][3] = 0;
    blosum62[2][4] = -2;
    blosum62[2][5] = -3;
    blosum62[2][6] = 1;
    blosum62[2][7] = 0;
    blosum62[2][8] = -2;
    blosum62[2][9] = 0;
    blosum62[2][10] = -3;
    blosum62[2][11] = -2;
    blosum62[2][12] = 2;
    blosum62[2][13] = -1;
    blosum62[2][14] = -3;
    blosum62[2][15] = -2;
    blosum62[2][16] = -1;
    blosum62[2][17] = -1;
    blosum62[2][18] = -3;
    blosum62[2][19] = -2;
    blosum62[2][20] = -3;
    blosum62[2][21] = -1;
    blosum62[2][22] = 0;
    blosum62[2][23] = -1;
    blosum62[3][1] = -2;
    blosum62[3][2] = 0;
    blosum62[3][3] = 6;
    blosum62[3][4] = 1;
    blosum62[3][5] = -3;
    blosum62[3][6] = 0;
    blosum62[3][7] = 0;
    blosum62[3][8] = 0;
    blosum62[3][9] = 1;
    blosum62[3][10] = -3;
    blosum62[3][11] = -3;
    blosum62[3][12] = 0;
    blosum62[3][13] = -2;
    blosum62[3][14] = -3;
    blosum62[3][15] = -2;
    blosum62[3][16] = 1;
    blosum62[3][17] = 0;
    blosum62[3][18] = -4;
    blosum62[3][19] = -2;
    blosum62[3][20] = -3;
    blosum62[3][21] = 3;
    blosum62[3][22] = 0;
    blosum62[3][23] = -1;
    blosum62[4][1] = -2;
    blosum62[4][2] = -2;
    blosum62[4][3] = 1;
    blosum62[4][4] = 6;
    blosum62[4][5] = -3;
    blosum62[4][6] = 0;
    blosum62[4][7] = 2;
    blosum62[4][8] = -1;
    blosum62[4][9] = -1;
    blosum62[4][10] = -3;
    blosum62[4][11] = -4;
    blosum62[4][12] = -1;
    blosum62[4][13] = -3;
    blosum62[4][14] = -3;
    blosum62[4][15] = -1;
    blosum62[4][16] = 0;
    blosum62[4][17] = -1;
    blosum62[4][18] = -4;
    blosum62[4][19] = -3;
    blosum62[4][20] = -3;
    blosum62[4][21] = 4;
    blosum62[4][22] = 1;
    blosum62[4][23] = -1;
    blosum62[5][1] = 0;
    blosum62[5][2] = -3;
    blosum62[5][3] = -3;
    blosum62[5][4] = -3;
    blosum62[5][5] = 9;
    blosum62[5][6] = -3;
    blosum62[5][7] = -4;
    blosum62[5][8] = -3;
    blosum62[5][9] = -3;
    blosum62[5][10] = -1;
    blosum62[5][11] = -1;
    blosum62[5][12] = -3;
    blosum62[5][13] = -1;
    blosum62[5][14] = -2;
    blosum62[5][15] = -3;
    blosum62[5][16] = -1;
    blosum62[5][17] = -1;
    blosum62[5][18] = -2;
    blosum62[5][19] = -2;
    blosum62[5][20] = -1;
    blosum62[5][21] = -3;
    blosum62[5][22] = -3;
    blosum62[5][23] = -2;
    blosum62[6][1] = -1;
    blosum62[6][2] = 1;
    blosum62[6][3] = 0;
    blosum62[6][4] = 0;
    blosum62[6][5] = -3;
    blosum62[6][6] = 5;
    blosum62[6][7] = 2;
    blosum62[6][8] = -2;
    blosum62[6][9] = 0;
    blosum62[6][10] = -3;
    blosum62[6][11] = -2;
    blosum62[6][12] = 1;
    blosum62[6][13] = 0;
    blosum62[6][14] = -3;
    blosum62[6][15] = -1;
    blosum62[6][16] = 0;
    blosum62[6][17] = -1;
    blosum62[6][18] = -2;
    blosum62[6][19] = -1;
    blosum62[6][20] = -2;
    blosum62[6][21] = 0;
    blosum62[6][22] = 3;
    blosum62[6][23] = -1;
    blosum62[7][1] = -1;
    blosum62[7][2] = 0;
    blosum62[7][3] = 0;
    blosum62[7][4] = 2;
    blosum62[7][5] = -4;
    blosum62[7][6] = 2;
    blosum62[7][7] = 5;
    blosum62[7][8] = -2;
    blosum62[7][9] = 0;
    blosum62[7][10] = -3;
    blosum62[7][11] = -3;
    blosum62[7][12] = 1;
    blosum62[7][13] = -2;
    blosum62[7][14] = -3;
    blosum62[7][15] = -1;
    blosum62[7][16] = 0;
    blosum62[7][17] = -1;
    blosum62[7][18] = -3;
    blosum62[7][19] = -2;
    blosum62[7][20] = -2;
    blosum62[7][21] = 1;
    blosum62[7][22] = 4;
    blosum62[7][23] = -1;
    blosum62[8][1] = 0;
    blosum62[8][2] = -2;
    blosum62[8][3] = 0;
    blosum62[8][4] = -1;
    blosum62[8][5] = -3;
    blosum62[8][6] = -2;
    blosum62[8][7] = -2;
    blosum62[8][8] = 6;
    blosum62[8][9] = -2;
    blosum62[8][10] = -4;
    blosum62[8][11] = -4;
    blosum62[8][12] = -2;
    blosum62[8][13] = -3;
    blosum62[8][14] = -3;
    blosum62[8][15] = -2;
    blosum62[8][16] = 0;
    blosum62[8][17] = -2;
    blosum62[8][18] = -2;
    blosum62[8][19] = -3;
    blosum62[8][20] = -3;
    blosum62[8][21] = -1;
    blosum62[8][22] = -2;
    blosum62[8][23] = -1;
    blosum62[9][1] = -2;
    blosum62[9][2] = 0;
    blosum62[9][3] = 1;
    blosum62[9][4] = -1;
    blosum62[9][5] = -3;
    blosum62[9][6] = 0;
    blosum62[9][7] = 0;
    blosum62[9][8] = -2;
    blosum62[9][9] = 8;
    blosum62[9][10] = -3;
    blosum62[9][11] = -3;
    blosum62[9][12] = -1;
    blosum62[9][13] = -2;
    blosum62[9][14] = -1;
    blosum62[9][15] = -2;
    blosum62[9][16] = -1;
    blosum62[9][17] = -2;
    blosum62[9][18] = -2;
    blosum62[9][19] = 2;
    blosum62[9][20] = -3;
    blosum62[9][21] = 0;
    blosum62[9][22] = 0;
    blosum62[9][23] = -1;
    blosum62[10][1] = -1;
    blosum62[10][2] = -3;
    blosum62[10][3] = -3;
    blosum62[10][4] = -3;
    blosum62[10][5] = -1;
    blosum62[10][6] = -3;
    blosum62[10][7] = -3;
    blosum62[10][8] = -4;
    blosum62[10][9] = -3;
    blosum62[10][10] = 4;
    blosum62[10][11] = 2;
    blosum62[10][12] = -3;
    blosum62[10][13] = 1;
    blosum62[10][14] = 0;
    blosum62[10][15] = -3;
    blosum62[10][16] = -2;
    blosum62[10][17] = -1;
    blosum62[10][18] = -3;
    blosum62[10][19] = -1;
    blosum62[10][20] = 3;
    blosum62[10][21] = -3;
    blosum62[10][22] = -3;
    blosum62[10][23] = -1;
    blosum62[11][1] = -1;
    blosum62[11][2] = -2;
    blosum62[11][3] = -3;
    blosum62[11][4] = -4;
    blosum62[11][5] = -1;
    blosum62[11][6] = -2;
    blosum62[11][7] = -3;
    blosum62[11][8] = -4;
    blosum62[11][9] = -3;
    blosum62[11][10] = 2;
    blosum62[11][11] = 4;
    blosum62[11][12] = -2;
    blosum62[11][13] = 2;
    blosum62[11][14] = 0;
    blosum62[11][15] = -3;
    blosum62[11][16] = -2;
    blosum62[11][17] = -1;
    blosum62[11][18] = -2;
    blosum62[11][19] = -1;
    blosum62[11][20] = 1;
    blosum62[11][21] = -4;
    blosum62[11][22] = -3;
    blosum62[11][23] = -1;
    blosum62[12][1] = -1;
    blosum62[12][2] = 2;
    blosum62[12][3] = 0;
    blosum62[12][4] = -1;
    blosum62[12][5] = -3;
    blosum62[12][6] = 1;
    blosum62[12][7] = 1;
    blosum62[12][8] = -2;
    blosum62[12][9] = -1;
    blosum62[12][10] = -3;
    blosum62[12][11] = -2;
    blosum62[12][12] = 5;
    blosum62[12][13] = -1;
    blosum62[12][14] = -3;
    blosum62[12][15] = -1;
    blosum62[12][16] = 0;
    blosum62[12][17] = -1;
    blosum62[12][18] = -3;
    blosum62[12][19] = -2;
    blosum62[12][20] = -2;
    blosum62[12][21] = 0;
    blosum62[12][22] = 1;
    blosum62[12][23] = -1;
    blosum62[13][1] = -1;
    blosum62[13][2] = -1;
    blosum62[13][3] = -2;
    blosum62[13][4] = -3;
    blosum62[13][5] = -1;
    blosum62[13][6] = 0;
    blosum62[13][7] = -2;
    blosum62[13][8] = -3;
    blosum62[13][9] = -2;
    blosum62[13][10] = 1;
    blosum62[13][11] = 2;
    blosum62[13][12] = -1;
    blosum62[13][13] = 5;
    blosum62[13][14] = 0;
    blosum62[13][15] = -2;
    blosum62[13][16] = -1;
    blosum62[13][17] = -1;
    blosum62[13][18] = -1;
    blosum62[13][19] = -1;
    blosum62[13][20] = 1;
    blosum62[13][21] = -3;
    blosum62[13][22] = -1;
    blosum62[13][23] = -1;
    blosum62[14][1] = -2;
    blosum62[14][2] = -3;
    blosum62[14][3] = -3;
    blosum62[14][4] = -3;
    blosum62[14][5] = -2;
    blosum62[14][6] = -3;
    blosum62[14][7] = -3;
    blosum62[14][8] = -3;
    blosum62[14][9] = -1;
    blosum62[14][10] = 0;
    blosum62[14][11] = 0;
    blosum62[14][12] = -3;
    blosum62[14][13] = 0;
    blosum62[14][14] = 6;
    blosum62[14][15] = -4;
    blosum62[14][16] = -2;
    blosum62[14][17] = -2;
    blosum62[14][18] = 1;
    blosum62[14][19] = 3;
    blosum62[14][20] = -1;
    blosum62[14][21] = -3;
    blosum62[14][22] = -3;
    blosum62[14][23] = -1;
    blosum62[15][1] = -1;
    blosum62[15][2] = -2;
    blosum62[15][3] = -2;
    blosum62[15][4] = -1;
    blosum62[15][5] = -3;
    blosum62[15][6] = -1;
    blosum62[15][7] = -1;
    blosum62[15][8] = -2;
    blosum62[15][9] = -2;
    blosum62[15][10] = -3;
    blosum62[15][11] = -3;
    blosum62[15][12] = -1;
    blosum62[15][13] = -2;
    blosum62[15][14] = -4;
    blosum62[15][15] = 7;
    blosum62[15][16] = -1;
    blosum62[15][17] = -1;
    blosum62[15][18] = -4;
    blosum62[15][19] = -3;
    blosum62[15][20] = -2;
    blosum62[15][21] = -2;
    blosum62[15][22] = -1;
    blosum62[15][23] = -2;
    blosum62[16][1] = 1;
    blosum62[16][2] = -1;
    blosum62[16][3] = 1;
    blosum62[16][4] = 0;
    blosum62[16][5] = -1;
    blosum62[16][6] = 0;
    blosum62[16][7] = 0;
    blosum62[16][8] = 0;
    blosum62[16][9] = -1;
    blosum62[16][10] = -2;
    blosum62[16][11] = -2;
    blosum62[16][12] = 0;
    blosum62[16][13] = -1;
    blosum62[16][14] = -2;
    blosum62[16][15] = -1;
    blosum62[16][16] = 4;
    blosum62[16][17] = 1;
    blosum62[16][18] = -3;
    blosum62[16][19] = -2;
    blosum62[16][20] = -2;
    blosum62[16][21] = 0;
    blosum62[16][22] = 0;
    blosum62[16][23] = 0;
    blosum62[17][1] = 0;
    blosum62[17][2] = -1;
    blosum62[17][3] = 0;
    blosum62[17][4] = -1;
    blosum62[17][5] = -1;
    blosum62[17][6] = -1;
    blosum62[17][7] = -1;
    blosum62[17][8] = -2;
    blosum62[17][9] = -2;
    blosum62[17][10] = -1;
    blosum62[17][11] = -1;
    blosum62[17][12] = -1;
    blosum62[17][13] = -1;
    blosum62[17][14] = -2;
    blosum62[17][15] = -1;
    blosum62[17][16] = 1;
    blosum62[17][17] = 5;
    blosum62[17][18] = -2;
    blosum62[17][19] = -2;
    blosum62[17][20] = 0;
    blosum62[17][21] = -1;
    blosum62[17][22] = -1;
    blosum62[17][23] = 0;
    blosum62[18][1] = -3;
    blosum62[18][2] = -3;
    blosum62[18][3] = -4;
    blosum62[18][4] = -4;
    blosum62[18][5] = -2;
    blosum62[18][6] = -2;
    blosum62[18][7] = -3;
    blosum62[18][8] = -2;
    blosum62[18][9] = -2;
    blosum62[18][10] = -3;
    blosum62[18][11] = -2;
    blosum62[18][12] = -3;
    blosum62[18][13] = -1;
    blosum62[18][14] = 1;
    blosum62[18][15] = -4;
    blosum62[18][16] = -3;
    blosum62[18][17] = -2;
    blosum62[18][18] = 11;
    blosum62[18][19] = 2;
    blosum62[18][20] = -3;
    blosum62[18][21] = -4;
    blosum62[18][22] = -3;
    blosum62[18][23] = -2;
    blosum62[19][1] = -2;
    blosum62[19][2] = -2;
    blosum62[19][3] = -2;
    blosum62[19][4] = -3;
    blosum62[19][5] = -2;
    blosum62[19][6] = -1;
    blosum62[19][7] = -2;
    blosum62[19][8] = -3;
    blosum62[19][9] = 2;
    blosum62[19][10] = -1;
    blosum62[19][11] = -1;
    blosum62[19][12] = -2;
    blosum62[19][13] = -1;
    blosum62[19][14] = 3;
    blosum62[19][15] = -3;
    blosum62[19][16] = -2;
    blosum62[19][17] = -2;
    blosum62[19][18] = 2;
    blosum62[19][19] = 7;
    blosum62[19][20] = -1;
    blosum62[19][21] = -3;
    blosum62[19][22] = -2;
    blosum62[19][23] = -1;
    blosum62[20][1] = 0;
    blosum62[20][2] = -3;
    blosum62[20][3] = -3;
    blosum62[20][4] = -3;
    blosum62[20][5] = -1;
    blosum62[20][6] = -2;
    blosum62[20][7] = -2;
    blosum62[20][8] = -3;
    blosum62[20][9] = -3;
    blosum62[20][10] = 3;
    blosum62[20][11] = 1;
    blosum62[20][12] = -2;
    blosum62[20][13] = 1;
    blosum62[20][14] = -1;
    blosum62[20][15] = -2;
    blosum62[20][16] = -2;
    blosum62[20][17] = 0;
    blosum62[20][18] = -3;
    blosum62[20][19] = -1;
    blosum62[20][20] = 4;
    blosum62[20][21] = -3;
    blosum62[20][22] = -2;
    blosum62[20][23] = -1;
    blosum62[21][1] = -2;
    blosum62[21][2] = -1;
    blosum62[21][3] = 3;
    blosum62[21][4] = 4;
    blosum62[21][5] = -3;
    blosum62[21][6] = 0;
    blosum62[21][7] = 1;
    blosum62[21][8] = -1;
    blosum62[21][9] = 0;
    blosum62[21][10] = -3;
    blosum62[21][11] = -4;
    blosum62[21][12] = 0;
    blosum62[21][13] = -3;
    blosum62[21][14] = -3;
    blosum62[21][15] = -2;
    blosum62[21][16] = 0;
    blosum62[21][17] = -1;
    blosum62[21][18] = -4;
    blosum62[21][19] = -3;
    blosum62[21][20] = -3;
    blosum62[21][21] = 4;
    blosum62[21][22] = 1;
    blosum62[21][23] = -1;
    blosum62[22][1] = -1;
    blosum62[22][2] = 0;
    blosum62[22][3] = 0;
    blosum62[22][4] = 1;
    blosum62[22][5] = -3;
    blosum62[22][6] = 3;
    blosum62[22][7] = 4;
    blosum62[22][8] = -2;
    blosum62[22][9] = 0;
    blosum62[22][10] = -3;
    blosum62[22][11] = -3;
    blosum62[22][12] = 1;
    blosum62[22][13] = -1;
    blosum62[22][14] = -3;
    blosum62[22][15] = -1;
    blosum62[22][16] = 0;
    blosum62[22][17] = -1;
    blosum62[22][18] = -3;
    blosum62[22][19] = -2;
    blosum62[22][20] = -2;
    blosum62[22][21] = 1;
    blosum62[22][22] = 4;
    blosum62[22][23] = -1;
    blosum62[23][1] = 0;
    blosum62[23][2] = -1;
    blosum62[23][3] = -1;
    blosum62[23][4] = -1;
    blosum62[23][5] = -2;
    blosum62[23][6] = -1;
    blosum62[23][7] = -1;
    blosum62[23][8] = -1;
    blosum62[23][9] = -1;
    blosum62[23][10] = -1;
    blosum62[23][11] = -1;
    blosum62[23][12] = -1;
    blosum62[23][13] = -1;
    blosum62[23][14] = -1;
    blosum62[23][15] = -2;
    blosum62[23][16] = 0;
    blosum62[23][17] = 0;
    blosum62[23][18] = -2;
    blosum62[23][19] = -1;
    blosum62[23][20] = -1;
    blosum62[23][21] = -1;
    blosum62[23][22] = -1;
    blosum62[23][23] = -1;
}