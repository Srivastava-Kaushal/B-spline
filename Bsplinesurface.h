#ifndef BSPLINESURFACE_H
#define BSPLINESURFACE_H

#include <Eigen/Dense>
#include <vector>
#include "basisfunction.h"
#include "point.h"

using namespace Eigen;

class Bsplinesurface {
private:
    VectorXd knotvector1;
    VectorXd knotvector2;
    std::vector<std::vector<point>> controlPoint;
    int k;
    int l;
    int n;
    int m;
public:

    Bsplinesurface() : k(0), l(0), n(0), m(0) {
        // Initialize other member variables or leave them empty as needed.
    }

    Bsplinesurface(const VectorXd& knots1, const VectorXd& knots2, const std::vector<std::vector<point>>& controlPts, int k, int l, int n, int m)
        : knotvector1(knots1), knotvector2(knots2), controlPoint(controlPts), k(k), l(l), n(n), m(m) {}


    point PointCalculation(double u, double w) {
        point result{ 0.0f, 0.0f, 0.0f };

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                double basis1 = basisfunction(u, i, k, knotvector1);
                double basis2 = basisfunction(w, j, l, knotvector2);
                result.x += controlPoint[i][j].x * basis1 * basis2;
                result.y += controlPoint[i][j].y * basis1 * basis2;
                result.z += controlPoint[i][j].z * basis1 * basis2;
            }
        }
        return result;
    }
};
#endif // BSPLINE_H
