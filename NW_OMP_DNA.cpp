#include<iostream>
#include<fstream>
#include<algorithm>
#include<cmath>
#include<string>
#include<vector>
#include<omp.h>
#include<time.h>
//#include<iomanip>

using namespace std;

class myMatrix {
private:
string A;
string B;
int lenA;
int lenB;
int threadNums;
int matchScore;
int mismatchPenalty;
int gapPenalty;
vector<vector<int>> matrix;

public:
myMatrix(string &x, string &y, int thread, int match, int mismatch, int gap) {
    A = x;
    B = y;
    lenA = A.size();
    lenB = B.size();
    threadNums = thread;
    matchScore = match;
    mismatchPenalty = mismatch;
    gapPenalty = gap;
    for (int i = 0; i <= lenA; i++) {
        matrix.emplace_back(vector<int>(lenB+1,0));
    }
}
void initialize() {
    omp_set_num_threads(threadNums);
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i <= lenA; i++) {
            matrix[i][0] = i * gapPenalty;
        }
        #pragma omp for
        for (int j = 0; j <= lenB; j++) {
            matrix[0][j] = j * gapPenalty;
        }
    }
}
void fillElement(int i, int j) {
    if (i > lenA || j > lenB || i == 0 || j == 0) return; //comfirm that i and j are valid indices
    int match = (A[i-1] == B[j-1]) ? (matrix[i-1][j-1] + matchScore) : (matrix[i-1][j-1] + mismatchPenalty);
    int del = matrix[i-1][j] + gapPenalty;
    int ins = matrix[i][j-1] + gapPenalty;
    matrix[i][j] = max(match, max(del, ins));
}
void fillMatrix() {
    int cycle = ceil((float)lenA/threadNums); //separation of the matrix
    for (int k = 0; k < cycle; k++) {
        int start = k * threadNums + 1;
        int thread = (k == cycle - 1) ? (lenA - start + 1) : threadNums;
        int end = start + thread - 1;
        //sequential part 1
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
        //parallel part
        for (int count = 0; count < lenB - (thread - 1); count++) {
            omp_set_num_threads(thread);
            #pragma omp parallel for
            for (int i = start; i <= end; i++) {
                fillElement(i, thread - (i - start) + count);
            }
        }
        //sequential part 2
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
void backtracking(string &alignmentA, string &alignmentB) {
    int i = lenA;
    int j = lenB;
    while(i > 0 && j > 0) {
        int score = matrix[i][j];
        int scoreDiag = matrix[i-1][j-1];
        int scoreUp = matrix[i][j-1];
        int scoreLeft = matrix[i-1][j];
        int S = (A[i-1] == B[j-1]) ? matchScore : mismatchPenalty;
        if (score == scoreDiag + S) {
            alignmentA = A[i-1] + alignmentA;
            alignmentB = B[j-1] + alignmentB;
            i--;
            j--;
        }
        else if (score == scoreLeft + gapPenalty) {
            alignmentA = A[i-1] + alignmentA;
            alignmentB = "-" + alignmentB;
            i--;
        }
        else {
            alignmentA = "-" + alignmentA;
            alignmentB = B[j-1] + alignmentB;
            j--;
        }
    }
    while(i > 0) {
        alignmentA = A[i-1] + alignmentA;
        alignmentB = "-" + alignmentB;
        i--;
    }
    while(j > 0) {
        alignmentA = "-" + alignmentA;
        alignmentB = B[j-1] + alignmentB;
        j--;
    }
}
};

vector<string> readSequences(int lenA, int lenB) {
    ifstream f("data/5K_Sequence.fasta");
    int cycleA = lenA / 200; //the data length in one line is 200 in this file
    int cycleB = lenB / 200;
    string temp;
    char space;
    vector<string> ans(2, string(0, ' '));
    getline(f, temp);
    for (int i = 0; i < cycleA; i++) {
        getline(f, temp);
        getline(f, temp);
        f >> space; //remove the blank space at the beginning of each data line
        getline(f, temp);
        ans[0].append(temp);
    }
    for (int i = 0; i < cycleB; i++) {
        getline(f, temp);
        getline(f, temp);
        f >> space;
        getline(f, temp);
        ans[1].append(temp);
    }
    f.close();
    return ans;
}

int main() {
int lenA;
int lenB;
cout << "length of sequence A: ";
cin >> lenA;
cout << "length of sequence B: ";
cin >> lenB;
cout << "\nThe sequence data is being read from the file..." << endl;
vector<string> sequences = readSequences(lenA, lenB);
cout << "Reading has finished." << endl << endl;
int threads;
cout << "number of threads: ";
cin >> threads;
int match = 1;
int mismatch = -1;
int gap = -2;
int isDefault = 1;
cout << "Use default Score&Penalty value? (0: manually set; 1: default): ";
cin >> isDefault;
if (isDefault == 0) {
    cout << "match score: ";
    cin >> match;
    cout << "mismatch penalty: ";
    cin >> mismatch;
    cout << "gap penalty: ";
    cin >> gap;
}
cout << "\nstart to time..." << endl;
clock_t start = clock();
myMatrix dp(sequences[0], sequences[1], threads, match, mismatch, gap);
dp.initialize();
dp.fillMatrix();
string alignmentA;
string alignmentB;
dp.backtracking(alignmentA, alignmentB);
clock_t end = clock();
cout << "Execution time: " << (end - start) << "ms" << endl;
system("pause");
}