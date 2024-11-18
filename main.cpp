#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <Eigen/Dense>
// #include "eigen/Eigen/Dense"

#include "basisfunction.h"
#include "knotvector.h"
#include "chordlength.h"
#include "point.h"
#include "surfaceplot.h"
#include <functional>


using namespace std;
using namespace Eigen;

// Define parameters
const int POPULATION_SIZE = 10;
const int NUM_VARIABLES = 10;
const int MAX_GENERATIONS = 500;
const double MUTATION_RATE = 0.1;
const double MIN_VALUE = -5.0;
const double MAX_VALUE = 5.0;

// Function to generate random real-valued genomes
vector<double> generateRandomGenome() {
    vector<double> genome(NUM_VARIABLES);
    for (int i = 0; i < NUM_VARIABLES; ++i) {
        genome[i] = MIN_VALUE + static_cast<double>(rand()) / RAND_MAX * (MAX_VALUE - MIN_VALUE);
    }
    return genome;
}

// Function to calculate fitness of a genome (example: sphere function)
double calculateFitness(const vector<double>& genome) {
    double fitness = 0.0;
    for (double gene : genome) {
        fitness += gene * gene;
    }
    return fitness;
}

// Function to perform single-point crossover
vector<double> crossover(const vector<double>& parent1, const vector<double>& parent2) {
    vector<double> child(NUM_VARIABLES);
    int crossoverPoint = rand() % NUM_VARIABLES;
    for (int i = 0; i < crossoverPoint; ++i) {
        child[i] = parent1[i];
    }
    for (int i = crossoverPoint; i < NUM_VARIABLES; ++i) {
        child[i] = parent2[i];
    }
    return child;
}

// Function to perform mutation
void mutate(vector<double>& genome) {
    for (double& gene : genome) {
        if (rand() / static_cast<double>(RAND_MAX) < MUTATION_RATE) {
            gene += (static_cast<double>(rand()) / RAND_MAX - 0.5) * (MAX_VALUE - MIN_VALUE) * 0.1;
            // Mutation: add a small random value (-0.1 to 0.1) to the gene
            // This could be adjusted based on the problem domain
            if (gene < MIN_VALUE) gene = MIN_VALUE;
            if (gene > MAX_VALUE) gene = MAX_VALUE;
        }
    }
}


int main(int argc, char** argv) {
    // Get user input for order and number of control points
    int k, l, n, m;
    cout << "Enter Order in u-direction of B-Spline: ";
    cin >> k;
    cout << "Enter Order in v-direction of B-Spline: ";
    cin >> l;
    cout << "Enter Number of Control Points in u-direction (n): ";
    cin >> n;
    cout << "Enter Number of Control Points in v-direction (m): ";
    cin >> m;

    // Open the data file
    ifstream inputFile("data.txt");
    if (!inputFile.is_open()) {
        cerr << "Error: Unable to open data file." << endl;
        return 1;
    }

    // Read r and s from the data file
    int r, s;
    inputFile >> r >> s;
    MatrixXd datapoints(r * s, 3); // Matrix to store data points
    for (int i = 0; i < r * s; ++i) {
        inputFile >> datapoints(i, 0) >> datapoints(i, 1) >> datapoints(i, 2);
    }

    // Close the data file
    inputFile.close();

    // Displaying data points
    cout << "Data points matrix (r = " << r << ", s = " << s << "):" << endl << datapoints << endl;

    // Initialize knot vectors
    VectorXd knotvector_u = knotvector(n, k);
    cout << "Knot vector in u-direction: \n" << knotvector_u << endl;
    VectorXd knotvector_v = knotvector(m, l);
    cout << "Knot vector in v-direction: \n" << knotvector_v << endl;

    // Calculate parameters using chord length method
    VectorXd parameter_u = chordlength(datapoints, r, knotvector_u.maxCoeff());
    cout << "Parameter_u: \n" << parameter_u << endl;

    VectorXd parameter_v = chordlength(datapoints, s, knotvector_v.maxCoeff());
    cout << "Parameter_v: \n" << parameter_v << endl;

    // PSO Algorithm parameters
    srand(time(nullptr));

    // Initialize population
    vector<vector<double>> population(POPULATION_SIZE);
    for (int i = 0; i < POPULATION_SIZE; ++i) {
        population[i] = generateRandomGenome();
    }

    // Evolution loop
    int generation = 0;
    while (generation < MAX_GENERATIONS) {
        // Calculate fitness for each individual in the population
        vector<double> fitnessScores(POPULATION_SIZE);
        for (int i = 0; i < POPULATION_SIZE; ++i) {
            fitnessScores[i] = calculateFitness(population[i]);
        }

        // Find the best individual (lowest fitness)
        double bestFitness = fitnessScores[0];
        int bestIndex = 0;
        for (int i = 1; i < POPULATION_SIZE; ++i) {
            if (fitnessScores[i] < bestFitness) {
                bestFitness = fitnessScores[i];
                bestIndex = i;
            }
        }

        // Output the best fitness found in this generation
        cout << "In Generation #" << generation << " Best Fitness is " << bestFitness << endl;

        // Select parents for crossover (tournament selection)
        int parent1Index = rand() % POPULATION_SIZE;
        int parent2Index = rand() % POPULATION_SIZE;
        if (fitnessScores[parent1Index] < fitnessScores[parent2Index])
            parent2Index = parent1Index;

        // Perform crossover to create new offspring
        vector<double> child = crossover(population[parent1Index], population[parent2Index]);

        // Mutate the child
        mutate(child);

        // Replace the worst individual in the population with the child
        double worstFitness = fitnessScores[0];
        int worstIndex = 0;
        for (int i = 1; i < POPULATION_SIZE; ++i) {
            if (fitnessScores[i] > worstFitness) {
                worstFitness = fitnessScores[i];
                worstIndex = i;
            }
        }
        population[worstIndex] = child;

        generation++;
    }



    // Other parts of the code remain the same
    // Calculate matrices N and M, matrix C, and control points
    MatrixXd N(parameter_u.size(), n);
    for (int c = 0; c < parameter_u.size(); c++) {
        for (int i = 0; i < n; i++) {
            N(c, i) = basisfunction(parameter_u(c), i, k, knotvector_u);
        }
    }
    cout << "Matrix N: \n" << N << endl;

    MatrixXd M(parameter_v.size(), m);
    for (int d = 0; d < parameter_v.size(); d++) {
        for (int i = 0; i < m; i++) {
            M(d, i) = basisfunction(parameter_v(d), i, l, knotvector_v);
        }
    }
    cout << "Matrix M: \n" << M << endl;

    // Calculate matrix C
    int rows = r * s;
    int cols = n * m;
    MatrixXd C(rows, cols);
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < s; ++j) {
            for (int p = 0; p < n; ++p) {
                for (int q = 0; q < m; ++q) {
                    C(i * s + j, p * m + q) = N(i, p) * M(j, q);
                }
            }
        }
    }
    cout << "Matrix C: \n" << C << endl;

    // Calculate control points
    MatrixXd controlPoints;
    if (rows == cols) {
        controlPoints = C.colPivHouseholderQr().solve(datapoints);
    }
    else {
        controlPoints = (C.transpose() * C).colPivHouseholderQr().solve(C.transpose() * datapoints);
    }
    cout << "Control Points: \n" << controlPoints << endl;

    // Control points conversion
    vector<vector<point>> controlPoints1(n, vector<point>(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            point p1;
            p1.x = controlPoints(i * m + j, 0);
            p1.y = controlPoints(i * m + j, 1);
            p1.z = controlPoints(i * m + j, 2);
            controlPoints1[i][j] = p1;
        }
    }

    // Display control points
    cout << "\nDisplay of Control Points:" << endl;
    for (const auto& row : controlPoints1) {
        for (const auto& p : row) {
            cout << "(" << p.x << ", " << p.y << ", " << p.z << ") ";
        }
        cout << endl;
    }

    // Call graphics function (assuming it is available)
    graphics(argc, argv, datapoints, controlPoints1, knotvector_u, knotvector_v, k, l, n, m, r, s);

    return 0;
}
