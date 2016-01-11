This is an implementation of the nearest neighbour and 2-opt algorithm for approximating tours for instances of the TSP-problem.

## Nearest neighbour ##
The nearest neighbour algorithm is a constructive algorithm that procedurally builds a tour by visiting every point in the input data and creating an edge between the current edge and the edge not yet picked that forms the shortest path from the current point. The closest neighbour can be found by searching the row corresponding to the current point in the distance matrix for the smallest value. If this point is not yet in the tour, it is added, else the second shortest distance point is searched for.

## Greedy 2-opt ##
The 2-opt algorithm is an optimization algorithm that iterates through every edge in the tour, and for each edge i, it iterates through the rest of the edges and examines if the tour becomes shorter if the two edges i and j is replaced by two new edges created by the respective end points of the two old edges. This procedure is done over and over again until no further improvement can be made and the tour is so called 2-optimized.
