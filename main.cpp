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

    int totalDistance = 0;

    for (int i = 0; i < N - 1; i++) {
        totalDistance += euclideanDistance((*X)[i], (*Y)[i], (*X)[i + 1], (*Y)[i + 1]);
    }

    totalDistance += euclideanDistance((*X)[N - 1], (*Y)[N - 1], (*X)[0], (*Y)[0]);

    while (true) {
        change = false;
        for (int j = 0; j < N - 1; j++) {
            int prev = j - 1;
            if (prev < 0) {
                prev = N - 1;
            }
            int next = j + 1;
            int nextnext = j + 2;
            if (nextnext > N - 1) {
                nextnext = 0;
            }
            if (next > N - 1) {
                next = 0;
                nextnext = 1;
            }
            int newDistance = totalDistance
                - euclideanDistance((*X)[prev], (*Y)[prev], (*X)[j], (*Y)[j])
                - euclideanDistance((*X)[nextnext], (*Y)[nextnext], (*X)[next], (*Y)[next])
                + euclideanDistance((*X)[prev], (*Y)[prev], (*X)[next], (*Y)[next])
                + euclideanDistance((*X)[j], (*Y)[j], (*X)[nextnext], (*Y)[nextnext]);
            if (newDistance < totalDistance) {
                totalDistance = newDistance;
                swapEdges(path, j, next);
                swapEdges(X, j, next);
                swapEdges(Y, j, next);
                change = true;
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

    return 0;
}
