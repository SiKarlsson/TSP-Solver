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

void swapEdges(double* coordinates[][2], int a, int b)
{
    double temp[2];
    temp[0] = *coordinates[a][0];
    temp[1] = *coordinates[a][1];
    *coordinates[a] = *coordinates[b];
    *coordinates[b] = temp;
}

void printIndexes(double coordinates[][2], int N)
{
    for (int i = 0; i < N; ++i)
    {
        cout << i << endl;
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

    double* xlist = new double[N];

    string x;
    string y;

    for (int i = 0; i < N; ++i) {
        cin >> x >> y;

        (coordinates)[i][0] = stod(x);
        (coordinates)[i][1] = stod(y);

        if (printouts) {
            cout << "x: " << (coordinates)[i][0] << ", y: " << (coordinates)[i][1] << endl;
        }
    }

    printIndexes(coordinates, N);

    swapEdges(arp, 0, 1);

    printIndexes(coordinates, N);

    return 0;
}
