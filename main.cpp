#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <limits.h>
#include <algorithm>
#include <queue>

// toggle debug print outs
const bool DEBUG = false;

// type edge<A, B>
typedef std::pair<int, int> edge;

inline int euclideanDistance(double x1, double y1, double x2, double y2)
{
    double distSum = pow(x1 - x2, 2) + pow(y1 - y2, 2);

    return round(sqrt(distSum));
}

inline std::vector<int> swapEdges(std::vector<int> tour, int x, int y, std::vector<int> &pos)
{
    int N = tour.size();
    int numSwaps = (((x <= y ? y - x : (y + N) - x) + 1)/2);
    int i = x;
    int j = y;
    for (int n = 0; n < numSwaps; ++n) {
        std::swap(tour[i], tour[j]);
        pos[tour[i]] = i;
        pos[tour[j]] = j;
        i = (i + 1) % N;        // one step to middle from left
        j = ((j + N) - 1) % N;  // one step to middle from right
    }

    return tour;
}

inline std::vector<int> twoOpt(std::vector<int> tour,
                               std::vector<std::vector<int> > distances,
                               std::vector<std::vector<int> >neighbour,
                               std::vector<int> &pos,
                               int min, int max)
{
    int N = distances[0].size();    // Number of nodes
    int K = neighbour[0].size();    // Number of closest neighbours

    int a, b, c, d;                 // edges AB and CD
    int a_i, b_i, c_i, d_i;         // indexes of AB and CD

    bool twoOptimal = false;        // Instance is locally 2-optimal
    while (!twoOptimal) {
        twoOptimal = true;

        for (a_i = 0, b_i = 1; a_i < N; ++a_i, ++b_i) {
            a = tour[a_i];
            b = tour[b_i % N];

            for (int k = 0; k < K; ++k)
            {
                c_i = pos[neighbour[a][k]]; // position of k closest neighbour
                d_i = c_i + 1;              // position of k's next point
                c = tour[c_i];
                d = tour[d_i % N];

                // continue if looking at same pairs
                if (b == c || a == d) {
                    continue;
                }

                // threshold for swapping edges
                if (distances[a][c] + min > distances[a][b] + max) {
                    break;
                }

                // if swap reduces tour length
                if (
                    distances[a][c] + distances[b][d] <
                    distances[a][b] + distances[c][d]
                ) {
                    // swap edges between b_i and c_i
                    tour = swapEdges(tour, (b_i % N), c_i, pos);

                    // new max threshold
                    max = std::max(max, std::max(distances[a][c], distances[b][d]));

                    // 2-opt was found
                    twoOptimal = false;

                    break;
                }
            }
        }
    }

    return tour;
}

/*
inline std::vector<int> greedyTour(std::vector<int> tour, std::vector<std::vector<int> > distances)
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
*/

inline std::vector<int> greedyTour(std::vector<int> tour, std::vector<std::vector<int>> distances, std::vector<std::vector<int>> neighbour)
{
    std::vector<int> degrees(distances.size());
    for (int i = 0; i < degrees.size(); ++i)
    {
        degrees[i] = 0;
    }
    std::vector<edge> edges;
    int x;
    int y;
    int minDistance;
    int dist;
    std::cout << "1" << std::endl;
    for (int n = 0; n < tour.size(); ++n) {
        std::cout << "n: " << n << std::endl;
        minDistance = INT_MAX;
        for (int i = 0; i < distances.size(); ++i) {
            dist = distances[i][neighbour[i][0]];

            if (dist < minDistance && degrees[i] < 2 && degrees[neighbour[i][0]] < 2) {
                minDistance = dist;
                x = i;
                y = neighbour[i][0];
            }

        }
        std::cout << "x: " << x << ", y: " << y << std::endl;
        degrees[x]++;
        degrees[y]++;
        neighbour[x].erase(neighbour[x].begin());
        int pos = find(neighbour[y].begin(), neighbour[y].end(), x) - neighbour[y].begin();
        neighbour[y].erase(neighbour[y].begin() + pos);
        //edges.push_back(edge(std::min(x, y), std::max(x, y)));
        edges.push_back(edge(x, y));

        /*
        for (int i = 0; i < neighbour.size(); ++i) {
            for (int j = 0; j < neighbour[i].size(); ++j) {
                std::cout << neighbour[i][j] << " ";
            }
            std::cout << std::endl;
        }
        */
    }
    std::cout << "2" << std::endl;

    int test;

    tour[0] = edges[0].first;
    tour[1] = edges[0].second;
    test = edges[0].second;
    edges[0] = edge(-1, -1);

    std::cout << "3" << std::endl;
    for (int i = 2; i < edges.size(); ++i) {
        for (int j = 0; j < edges.size(); ++j)
        {
            if (edges[j].first == test) {
                tour[i] = edges[j].second;
                test = edges[j].second;
                edges[j] = edge(-1, -1);
                break;
            } else if (edges[j].second == test ) {
                tour[i] = edges[j].first;
                test = edges[j].first;
                edges[j] = edge(-1, -1);
                break;
            }
        }
    }
    std::cout << "4" << std::endl;

    return tour;
}

inline std::vector<int> nearestNeighbour(std::vector<int> tour, std::vector<std::vector<int>> distances)
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

int main()
{
    // number of points
    int N;
    std::cin >> N;

    // number of neighbours to check on 2opt
    int K = 80;
    K = std::min(K, N);

    int numIt = 1;

    // shortest distance between two points
    int min = INT_MAX;

    // longest distance between two points
    int max = 0;

    // tour of points visited
    std::vector<int> tour(N);

    // indexes of points
    std::vector<int> index(N);

    // initial points vectors for X and Y coordinates
    std::vector<double> X (N);
    std::vector<double> Y (N);

    // distance matrix
    std::vector<std::vector<int>> distances (N, std::vector<int> (N));

    // neighbour matrix
    std::vector<std::vector<int>> neighbour (N);

    for (int i = 0; i < N; ++i)
    {
        // i:th line of input contains X and Y coordinates of point i
        std::cin >> X[i] >> Y[i];
    }

    for (int i = 0; i < N; ++i)
    {
        // initial tour
        tour[i] = i;

        for (int j = (i + 1); j < N; ++j)
        {
            // eucledian distance of (i, j) = (j, i)
            distances[i][j] = distances[j][i] = euclideanDistance(
                X[i], Y[i],
                X[j], Y[j]
            );

            min = std::min(min, distances[i][j]);
            max = std::max(max, distances[i][j]);
        }

        index[tour[i]] = i;
    }

    // compare distances of two edges
    auto compare = [&](edge a, edge b) {
        return distances[a.first][a.second] > distances[b.first][b.second];
    };

    for (int i = 0; i < N; ++i)
    {
        std::vector<edge> edges(N);

        // fill edge vector
        for (int j = 0; j < N; ++j)
        {
            edges[j] = edge(i, j);
        }

        // add edges to priority queue in ascending order using an internal heap
        std::priority_queue<edge, std::vector<edge>, decltype(compare)> queue {compare, std::move(edges)};

        for (int j = 0; j < K; ++j)
        {
            // j closest edge to i
            edge test = queue.top();
            if (test.first != test.second) {
                neighbour[i].push_back(test.second);
                if (DEBUG) {
                    std::cout << "[" << test.first << ", " << test.second << "] = " << distances[test.first][test.second] << " ";
                }
            }
            // pop edge on queue
            queue.pop();
        }
        if (DEBUG) {
            std::cout << std::endl;
        }
    }

    std::vector<int> bestTour(N);
    int minDistance = INT_MAX;

    for (int i = 0; i < numIt; i++) {
        // initial nearest neighbour serach
        //tour = nearestNeighbour(tour, distances);

        tour = greedyTour(tour, distances, neighbour);

        // optimize with 2-opt limited to K neighbours
        //tour = twoOpt(tour, distances, neighbour, index, min, max);

        int totalDistance = 0;
        for (int i = 0; i < tour.size() - 1; ++i) {
            totalDistance += distances[tour[i]][tour[i+1]];
        }
        totalDistance += distances[tour[tour.size()-1]][tour[0]];

        if (totalDistance < minDistance) {
            minDistance = totalDistance;
            bestTour = tour;
        }
    }

    // print final tour
    for (int i = 0; i < bestTour.size(); ++i)
    {
        std::cout << bestTour[i] << std::endl;
    }

    if (DEBUG) {
        int totalDistance = 0;
        for (int i = 0; i < tour.size() - 1; ++i)
        {
            std::cout << X[bestTour[i]] << " " << Y[bestTour[i]] << std::endl;
            totalDistance += distances[tour[i]][tour[i+1]];
        }
        totalDistance += distances[tour[tour.size()-1]][tour[0]];

        std::cout << "total: " << totalDistance << std::endl;
    }

    return 0;
}
