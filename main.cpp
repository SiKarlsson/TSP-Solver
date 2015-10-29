#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <limits.h>

using namespace std;

const bool printouts = false;

int euclideanDistance(double x1, double y1, double x2, double y2)
{
    double distSum = pow(x1 - x2, 2) + pow(y1 - y2, 2);

    return round(sqrt(distSum));
}

void swapEdges(double *I[], int a, int b)
{
    double temp = (*I)[a];
    (*I)[a] = (*I)[b];
    (*I)[b] = temp;
}

void printIndexes(double I[], int N)
{
    for (int i = 0; i < N; ++i)
    {
        cout << I[i] << endl;
    }
}

void twoOpt(int N, double *X[], double *Y[], double *path[])
{
    bool change;

    while (true) {
        change = false;
        for (int j = 0; j < N - 1; ++j) {
            for (int k = (j + 1); k < N - 1; ++k) {
                //cout << "j: " << j << ", k: " << k << endl;
                int prev = j - 1;
                if (prev < 0) {
                    prev = N - 1;
                }
                int next = k + 1;
                if (next > N - 1) {
                    next = 0;
                }

                int new_edge = euclideanDistance((*X)[j], (*Y)[j], (*X)[next], (*Y)[next])
                   + euclideanDistance((*X)[k], (*Y)[k], (*X)[prev], (*Y)[prev]);

                //cout << "new_edge: " << new_edge << endl;

                int old_edge = euclideanDistance((*X)[prev], (*Y)[prev], (*X)[j], (*Y)[j])
                   + euclideanDistance((*X)[k], (*Y)[k], (*X)[next], (*Y)[prev]);

                //cout << "old_edge: " << old_edge << endl;

                if (new_edge
                  < old_edge) {
                    //cout << "change: true" << endl;
                    swapEdges(path, j, k);
                    swapEdges(X, j, k);
                    swapEdges(Y, j, k);
                    //change = true;
                    break;
                }
            }
            if (change) {
                break;
            }
        }
        if (change) {
            continue;
        } else {
            break;
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

    twoOpt(N, &X, &Y, &path);

    printIndexes(path, N);


    int totalDistance = 0;

    for (int i = 0; i < N - 1; i++) {
        totalDistance += euclideanDistance(X[i], Y[i], X[i + 1], Y[i + 1]);
    }

    totalDistance += euclideanDistance(X[N - 1], Y[N - 1], X[0], Y[0]);

    cout << "tot: " << totalDistance << endl;

    return 0;
}
