VectorXd chordlength(const MatrixXd& datapoints, int size, double range) {
    VectorXd u(size * size);
    u.setZero(); // Initialize u vector with zeros

    for (int j = 0; j < size; j++) {
        double s1 = 0, s2 = 0;
        for (int i = 1; i < size; i++) {
            double dx = datapoints(j * size + i, 0) - datapoints(j * size + i - 1, 0);
            double dy = datapoints(j * size + i, 1) - datapoints(j * size + i - 1, 1);
            double dz = datapoints(j * size + i, 2) - datapoints(j * size + i - 1, 2);
            u(j * size + i) = sqrt(dx * dx + dy * dy + dz * dz);
            s1 += u(j * size + i);
        }
        for (int ii = 0; ii < size; ii++) {
            s2 += u(j * size + ii);
            u(j * size + ii) = (s2 / s1) * range;
        }
    }

    // Adjusting the value slightly smaller than 1
    for (int i = 0; i < u.size(); ++i) {
        if (u(i) >= 1.0) {
            u(i) -= 0.000001;
        }
    }

    return u;
}
