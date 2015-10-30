#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <limits.h>
#include <algorithm>

using namespace std;

const bool DEBUG = false;

std::vector<int> tour;
std::vector<double> X;
std::vector<double> Y;
std::vector<std::vector<int> > distances;

inline int euclideanDistance(double x1, double y1, double x2, double y2)
{
    double distSum = pow(x1 - x2, 2) + pow(y1 - y2, 2);

    return round(sqrt(distSum));
}

inline std::vector<double> reverseEdges(std::vector<double> TEST, int a, int b)
{
    int temp;
    int count = abs(b-a)+1;
    int mini = min(a, b);
    for (int i = 0; i < count/2; ++i) {
        temp = (TEST)[mini+count-i-1];
        (TEST)[mini+count-i-1] = (TEST)[mini+i];
        (TEST)[mini+i] = temp;
    }

    return TEST;

}

inline std::vector<int> swapEdges(std::vector<int> TEST, int a, int b)
{
    int temp;
    int count = abs(b-a)+1;
    int mini = min(a, b);
    for (int i = 0; i < count/2; ++i) {
        temp = (TEST)[mini+count-i-1];
        (TEST)[mini+count-i-1] = (TEST)[mini+i];
        (TEST)[mini+i] = temp;
    }

    return TEST;
}

inline std::vector<int> twoOpt(int N)
{
    if (N < 4) {
        return tour;
    }

    bool change = true;

    int times = 200;

    int ed, mini, maxi, new_edge, old_edge;

    while (change && times > 0) {
    //while (change) {
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

                new_edge = (distances)[tour[j]][tour[next]] + (distances)[tour[k]][tour[prev]];
                old_edge = (distances)[tour[j]][tour[prev]] + (distances)[tour[k]][tour[next]];

                if (new_edge < old_edge) {
                    tour = swapEdges(tour, j, k);
                    /*
                    X = reverseEdges(X, j, k);
                    Y = reverseEdges(Y, j, k);

                    mini = min(j, k);
                    maxi = max(j, k);

                    for (int l = mini; l <= maxi; l++) {
                        for (int m = 0; m < N; m++) {
                            if (l != m) {
                                ed = euclideanDistance((X)[l], (Y)[l], (X)[m], (Y)[m]);
                                (distances)[l][m] = ed;
                                (distances)[m][l] = ed;
                            }
                        }
                    }
                    */

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

    return tour;
}

inline std::vector<double> swapDoubles(std::vector<double> v, double x, double y)
{
    int temp = v[x];
    v[x] = v[y];
    v[y] = temp;

    return v;
}

std::vector<int> nearestNeighbour()
{
    int nearestIndex, nearestFound, temp;

    for (int i = 0; i < (tour.size() - 1); ++i)
    {
        nearestFound = INT_MAX;
        nearestIndex = i + 1;

        for (int j = (i + 1); j < tour.size(); ++j)
        {
            if (distances[tour[i]][tour[j]] < nearestFound) {
                nearestIndex = j;
                nearestFound = distances[tour[i]][tour[j]];
            }
        }

        temp = tour[i + 1];
        tour[i + 1] = tour[nearestIndex];
        tour[nearestIndex] = temp;

        temp = X[i + 1];
        X[i + 1] = X[nearestIndex];
        X[nearestIndex] = temp;

        temp = Y[i + 1];
        Y[i + 1] = Y[nearestIndex];
        Y[nearestIndex] = temp;
    }

    return tour;
}

int main()
{
    int N;

    std::cin >> N;

    X.resize(N);
    Y.resize(N);
    tour.resize(N);

    double x;
    double y;

    for (int i = 0; i < N; ++i) {
        std::cin >> x >> y;

        X[i] = x;
        Y[i] = y;
        tour[i] = i;
    }

    int ed;
    distances.resize(N);
    for (int i = 0; i < N; ++i)
    {
        distances[i].resize(N);
    }
    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            if (i != j) {
                ed = euclideanDistance(X[i], Y[i], X[j], Y[j]);
                distances[i][j] = ed;
                distances[j][i] = ed;
            }
        }
    }

    tour = nearestNeighbour();

    tour = twoOpt(N);

    for (int i = 0; i < tour.size(); ++i)
    {
        std::cout << tour[i] << endl;
    }

    return 0;
}
