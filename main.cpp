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

int main()
{
    int N;

    // read N
    cin >> N;

    if (printouts) {
        // print N
        cout << "N = " << N << endl;
    }

    double coordinates[N][2];

    string x;
    string y;

    for (int i = 0; i < N; ++i) {
        cin >> x >> y;

        coordinates[i][0] = stod(x);
        coordinates[i][1] = stod(y);

        if (printouts) {
            cout << "x: " << coordinates[i][0] << ", y: " << coordinates[i][1] << endl;
        }
    }

    return 0;
}
