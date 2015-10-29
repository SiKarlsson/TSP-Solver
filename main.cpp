#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <limits.h>
#include <algorithm>

using namespace std;

const bool printouts = false;

int euclideanDistance(double x1, double y1, double x2, double y2)
{
    double distSum = pow(x1 - x2, 2) + pow(y1 - y2, 2);

    return round(sqrt(distSum));
}

void swapEdges(double *I[], int a, int b)
{
    int temp;
    int count = abs(b-a)+1;
    int mini = min(a, b);
    for (int i = 0; i < count/2; ++i) {
        temp = (*I)[mini+count-i-1];
        (*I)[mini+count-i-1] = (*I)[mini+i];
        (*I)[mini+i] = temp;
    }
}

void printIndexes(double I[], int N)
{
    for (int i = 0; i < N; ++i)
    {
        cout << I[i] << endl;
    }
}

void twoOpt(int N, double *X[], double *Y[], double *path[], vector<vector<int> > *distances)
{
    if (N < 4) {
        return;
    }

    bool change = true;

    int times = 350;

    while (change && times > 0) {
        --times;
        change = false;
        for (int j = 0; j < N; ++j) {
            for (int k = (j + 1); k < N - 1; ++k) {
                int prev = j - 1;
                if (prev < 0) {
                    prev = N - 1;
                }
                int next = k + 1;
                if (next > N - 1) {
                    next = 0;
                }

                if (j == 0 && k == N - 2) {
                    prev = 1;
                    next = N-2;
                }

                int new_edge = (*distances)[j][next] + (*distances)[k][prev];
                int old_edge = (*distances)[j][prev] + (*distances)[k][next];

                if (new_edge < old_edge) {
                    swapEdges(path, j, k);
                    swapEdges(X, j, k);
                    swapEdges(Y, j, k);

                    int mini = min(j, k);
                    int maxi = max(j, k);
                    int ed;
                    
                    for (int l = mini; l <= maxi; l++) {
                        for (int m = 0; m < N; m++) {
                            if (l != m) {
                                ed = euclideanDistance((*X)[l], (*Y)[l], (*X)[m], (*Y)[m]);
                                (*distances)[l][m] = ed;
                                (*distances)[m][l] = ed;
                            }
                        }
                    }

                    change = true;
                }
                if (change) {
                    break;
                }
            }
            if (change) {
                break;
            }
        }
    }
}

int main()
{
    int N;

    cin >> N;

    if (printouts) {
        cout << "N = " << N << endl;
    }

    double* X = new double[N];
    double* Y = new double[N];
    double* path = new double[N];

    string x;
    string y;

    for (int i = 0; i < N; ++i) {
        cin >> x >> y;

        X[i] = stod(x);
        Y[i] = stod(y);
        path[i] = i;

        if (printouts) {
            cout << "x: " << X[i] << ", y: " << Y[i] << endl;
        }
    }

    vector<vector<int> > distances(N);
    for (int i = 0; i < N; i++) {
        distances[i].resize(N);
        for (int j = 0; j < N; j++) {
            if (i != j) {
                distances[i][j] = euclideanDistance(X[i], Y[i], X[j], Y[j]);
            }
        }
    }

    twoOpt(N, &X, &Y, &path, &distances);

    printIndexes(path, N);


    int totalDistance = 0;

    for (int i = 0; i < N - 1; i++) {
        totalDistance += euclideanDistance(X[i], Y[i], X[i + 1], Y[i + 1]);
    }

    totalDistance += euclideanDistance(X[N - 1], Y[N - 1], X[0], Y[0]);

    return 0;
}
