#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <limits.h>
#include <algorithm>
#include <queue>

using namespace std;

const bool DEBUG = true;

inline int euclideanDistance(double x1, double y1, double x2, double y2)
{
    double distSum = pow(x1 - x2, 2) + pow(y1 - y2, 2);

    return round(sqrt(distSum));
}

inline std::vector<int> swapEdges(std::vector<int> tour, int x, int y)
{
    int temp;
    int count = abs(y - x) + 1;
    int mini = min(x, y);
    for (int i = 0; i < (count / 2); ++i) {
        temp = tour[mini + count - i - 1];
        tour[mini + count - i - 1] = tour[mini + i];
        tour[mini + i] = temp;
    }

    return tour;
}

inline std::vector<int> twoOpt(std::vector<double> X, std::vector<double> Y, std::vector<int> tour, std::vector<std::vector<int> > distances)
{
    int N = tour.size();

    if (N < 4) {
        return tour;
    }

    bool change = true;

    int new_edge, old_edge, prev, next;

    while (change) {
        for (int j = 0; j < N; ++j) {
            for (int k = (j + 1); k < N - 1; ++k) {
                prev = j - 1;
                if (prev < 0) {
                    prev = N - 1;
                }
                next = k + 1;
                if (next > N - 1) {
                    next = 0;
                }

                if (j == 0 && k == N - 2) {
                    prev = 1;
                    next = N-2;
                }

                new_edge = distances[tour[j]][tour[next]] + distances[tour[k]][tour[prev]];
                old_edge = distances[tour[j]][tour[prev]] + distances[tour[k]][tour[next]];

                if (new_edge < old_edge) {
                    tour = swapEdges(tour, j, k);
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

inline std::vector<int> greedyTour(std::vector<int> tour, std::vector<double> X, std::vector<double> Y, std::vector<std::vector<int> > distances)
{
    int best;
    bool used[tour.size()];
    for (int i = 0; i < tour.size(); ++i)
    {
        used[i] = false;
    }
    used[0] = true;
    for (int i = 1; i < tour.size(); ++i) {
        best = -1;
        for (int j = 0; j < tour.size(); ++j) {
            if (! used[j] && (best == -1 || distances[tour[i - 1]][j] < distances[tour[i - 1]][best])) {
                best = j;
            }
        }
        tour[i] = best;
        used[best] = true;
    }

    return tour;
}

inline std::vector<int> nearestNeighbour(std::vector<double> X, std::vector<double> Y, std::vector<int> tour, std::vector<std::vector<int> > distances)
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
    }

    return tour;
}

inline std::vector<int> twoOptNn(std::vector<int> tour, std::vector<std::vector<int> > distances, std::vector<std::vector<int> > nn, int min, int max)
{
    int N = tour.size();
    int K = nn[0].size();

    if (N < 4) {
        return tour;
    }

    bool change = true;

    int new_edge, old_edge, prev, next;

    while (change) {
        change = false;
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < K; ++j)
            {
                prev = (i-1 < 0) ? N-1 : i-1;
                next = nn[i][j] + 1;
                next = (next > N - 1) ? 0 : next;

                if (tour[i] == tour[prev]) {
                    continue;
                }

                if (tour[i] == tour[nn[i][j]] || tour[next] == tour[prev]) {
                    continue;
                }

                if (distances[tour[prev]][tour[nn[i][j]]] + min > distances[tour[prev]][tour[i]] + max) {
                    break;
                }

                if (DEBUG) {
                    cout << "kollar på " << i << "-" << prev << " och " << nn[i][j] << "-" << next << endl;
                    //cout << "me värden " << tour[i] << "-" << tour[prev] << " och " << tour[nn[i][j]] << "-" << tour[next] << endl;
                }

                old_edge = distances[tour[prev]][tour[i]] + distances[tour[nn[i][j]]][tour[next]];
                new_edge = distances[tour[i]][tour[next]] + distances[tour[nn[i][j]]][tour[prev]];

                if (new_edge < old_edge) {
                    tour = swapEdges(tour, i, nn[i][j]);
                    change = true;
                    max = std::max(max, std::max(distances[tour[prev]][tour[nn[i][j]]], distances[tour[i]][tour[next]]));
                    break;
                }
            }
        }
    }

    return tour;
}

int main()
{
    int N;

    std::cin >> N;

    std::vector<int> tour(N);
    std::vector<double> X(N);
    std::vector<double> Y(N);
    std::vector<std::vector<int> > distances(N);

    double x, y;

    for (int i = 0; i < N; ++i) {
        std::cin >> x >> y;

        X[i] = x;
        Y[i] = y;
        tour[i] = i;
    }

    int ed;
    for (int i = 0; i < N; ++i)
    {
        distances[i].resize(N);
    }
    int min = INT_MAX;
    int max = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            if (i != j) {
                ed = euclideanDistance(X[i], Y[i], X[j], Y[j]);
                distances[i][j] = ed;
                distances[j][i] = ed;
                min = std::min(min, ed);
                max = std::max(max, ed);
            }
        }
    }

    if (DEBUG) {
        cout << "MIN DISTANCE = " << min << endl;
        cout << "MAX DISTANCE = " << max << endl;
    }

    int K = 20;
    K = std::min(K, N)-1;

    typedef pair<pair<int, int>, int> adj;

    auto comp = [](adj a, adj b) {
        return a.second > b.second;
    };

    std::vector<std::vector<int> > nn(N);

    for (int i = 0; i < N; ++i)
    {
        std::vector<adj> adjs(N);

        for (int j = 0; j < N; ++j)
        {
            adjs[j] = adj(pair<int, int>(i, j), distances[i][j]);
        }

        std::priority_queue<adj, std::vector<adj>, decltype(comp) > queue {comp, std::move(adjs)};

        for (int j = 0; j < (K + 1); ++j)
        {
            adj test = queue.top();
            if (test.first.first != test.first.second) {
                //nn[i].push_back(test.first.second);
                nn[i].push_back(test.first.second);
                if (DEBUG) {
                    cout << "[" << test.first.first << ", " << test.first.second << "] = " << test.second << " ";
                }
            }
            queue.pop();
        }
        if (DEBUG) {
            cout << endl;
        }
    }

    if (DEBUG) {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < K; ++j)
            {
                cout << nn[i][j] << " ";
            }
            cout << endl;
        }
    }

    //tour = greedyTour(tour, X, Y, distances);

    tour = nearestNeighbour(X, Y, tour, distances);

    //tour = twoOpt(X, Y, tour, distances);

    tour = twoOptNn(tour, distances, nn, min, max);

    for (int i = 0; i < tour.size(); ++i)
    {
        std::cout << tour[i] << endl;
    }

    if (DEBUG) {
        int totalDistance = 0;
        for (int i = 0; i < tour.size() - 1; ++i)
        {
            totalDistance += distances[tour[i]][tour[i+1]];
        }
        totalDistance += distances[tour[tour.size()-1]][tour[0]];

        std::cout << "total: " << totalDistance << endl;
    }

    return 0;
}
