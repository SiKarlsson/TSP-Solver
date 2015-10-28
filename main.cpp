#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

const bool printouts = true;

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
