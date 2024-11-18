#ifndef BASISFUNCTION_H
#define BASISFUNCTION_H

#include <Eigen/Dense>

using namespace Eigen;

double basisfunction(double t, int i, int k, const VectorXd& xi) {
    // t = parameter
    // k = order of B-Spline
    // xi = knot vector

    if (k == 1) {
        // Check if t is within the interval [xi(i), xi(i+1))
        if (xi[i] <= t && t < xi[i + 1]) {
            return 1.0;
        }
        else {
            return 0.0;
        }
    }
    else {
        // Recursive calculation of B-spline basis function
        double term1 = 0.0, term2 = 0.0;

        // Calculate the first term
        if ((xi[i + k - 1] - xi[i]) != 0) {
            term1 = (t - xi[i]) * basisfunction(t, i, k - 1, xi) / (xi[i + k - 1] - xi[i]);
        }

        // Calculate the second term
        if ((xi[i + k] - xi[i + 1]) != 0) {
            term2 = (xi[i + k] - t) * basisfunction(t, i + 1, k - 1, xi) / (xi[i + k] - xi[i + 1]);
        }

        // Combine terms to get the B-spline basis function
        return term1 + term2;
    }
}

#endif // BASISFUNCTION_H
