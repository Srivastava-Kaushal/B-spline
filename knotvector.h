#ifndef KNOTVECTOR_H
#define KNOTVECTOR_H

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

VectorXd knotvector(int n, int k) {
    std::cout << "CHOICE:\n1. Uniform\n2. Open Uniform\n3. Non-Uniform\n4. Exit\n";
    int choice;
    std::cin >> choice;

    VectorXd knot;
    switch (choice) {
        // Uniform knot vector
    case 1:
        knot.resize(n + k);
        for (int i = 0; i < n + k; ++i)
            knot(i) = i;
        break;

        // Open uniform knot vector
    case 2:
        knot.resize(n + k);
        for (int i = 1; i <= n + k; ++i) {
            if (i >= 1 && i <= k) {
                knot(i - 1) = 0;
            }
            else if (i >= k + 1 && i <= n) {
                knot(i - 1) = i - k;
            }
            else if (i >= n + 1 && i <= n + k) {
                knot(i - 1) = n - k + 1;
            }
        }
        break;

        // Non-uniform knot vector
    case 3:
        std::cout << "Enter the non-uniform knot vector:\n";
        knot.resize(n + k);
        for (int i = 0; i < n + k; ++i) {
            double val;
            std::cin >> val;
            knot(i) = val;
        }
        break;

    case 4:
        std::cout << "Exiting...\n";
        break;

    default:
        std::cout << "Wrong Choice!\n";
        break;
    }

    return knot;
}

#endif // KNOTVECTOR_H
