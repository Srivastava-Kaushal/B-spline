#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <GL/glut.h>
#include "Bsplinesurface.h"
#include "point.h"
#include "knotvector.h"

using namespace Eigen;

// Global variables
Bsplinesurface surface;
std::vector<std::vector<point>> controlPoints;
MatrixXd datapoints;
int n, m, r, s;

void initialize() {
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    gluLookAt(30.0, 30.0, 30.0, // camera position
        0.0, 0.0, 0.0,   // target position
        0.0, 1.0, 0.0); // up vector

    // Draw B-spline surface
    glColor3f(0.0f, 0.0f, 1.0f); // blue color
    double uStep = 0.05;
    double wStep = 0.05;

    for (double u = 0.0; u < 1.0; u += uStep) {
        glBegin(GL_LINE_STRIP);
        for (double w = 0.0; w < 2.0; w += wStep) {
            point pt = surface.PointCalculation(u, w);
            glVertex3d(pt.x, pt.y, pt.z);
        }
        glEnd();
    }

    for (double w = 0.0; w < 2.0; w += wStep) {
        glBegin(GL_LINE_STRIP);
        for (double u = 0.0; u < 1.0; u += uStep) {
            point pt = surface.PointCalculation(u, w);
            glVertex3d(pt.x, pt.y, pt.z);
        }
        glEnd();
    }

    // Draw control points
    glColor3f(0, 0, 0); // black color
    for (int i = 0; i < m; i++) {
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j < n; j++) {
            glVertex3f(controlPoints[i][j].x, controlPoints[i][j].y, controlPoints[i][j].z);
        }
        glEnd();
    }

    for (int j = 0; j < m; j++) {
        glBegin(GL_LINE_STRIP);
        for (int i = 0; i < n; i++) {
            glVertex3f(controlPoints[i][j].x, controlPoints[i][j].y, controlPoints[i][j].z);
        }
        glEnd();
    }

    // Draw data points
    glColor3f(1, 0, 0);
    glPointSize(2);
    for (int i = 0; i < r * s; i++) {
        glBegin(GL_POINTS);

        glVertex3f(datapoints(i, 0), datapoints(i, 1), datapoints(i, 2));

        glEnd();
    }

    glFlush();
    glutSwapBuffers();
}

void reshape(int width, int height) {
    glViewport(0, 0, static_cast<GLsizei>(width), static_cast<GLsizei>(height));
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, static_cast<GLdouble>(width) / static_cast<GLdouble>(height), 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

void graphics(int argc, char** argv, MatrixXd datapoint, std::vector<std::vector<point>> controlPoint, VectorXd knotvector1, VectorXd knotvector2, int k, int l, int n, int m, int r, int s) {
    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(1200, 800);
    glutCreateWindow("B-spline Surface");

    // Initialize OpenGL
    initialize();

    // Set up the B-spline surface with the given parameters
    surface = Bsplinesurface(knotvector1, knotvector2, controlPoint, k, l, n, m);
    controlPoints = controlPoint; // Assign control points to global variable
    datapoints = datapoint;
    ::n = n; // Assign n to global variable
    ::m = m; // Assign m to global variable
    ::r = r; // Assign r to global variable
    ::s = s; // Assign s to global variable

    // Register GLUT callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMainLoop();
}
