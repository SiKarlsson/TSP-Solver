#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>

using namespace std;

const bool printouts = true;

int euclideanDistance(double p1[], double p2[])
{
    double distSum = pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2);

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

int main()
{
    int N;

    // read N
    cin >> N;

    if (printouts) {
        // print N
        cout << "N = " << N << endl;
    }

    double* X = new double[N];
    double* Y = new double[N];
    double* I = new double[N];

    string x;
    string y;

    for (int i = 0; i < N; ++i) {
        cin >> x >> y;

        X[i] = stod(x);
        Y[i] = stod(y);
        I[i] = i;

        if (printouts) {
            cout << "x: " << X[i] << ", y: " << Y[i] << endl;
        }
    }

    /*
    for (int i = 1; i < N; ++i)
    {
        double distSum = pow(X[0] - p2[0], 2) + pow(p1[1] - p2[1], 2);

        return round(sqrt(distSum));
    }*/

    printIndexes(I, N);

    swapEdges(&I, 2, 1);

    printIndexes(I, N);

    return 0;
}
