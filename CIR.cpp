#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

Vector3d
KIR(Vector3d Wl_1n, Vector3d Wln, Vector3d Wl__1n, Matrix3d OmegaT, Matrix3d A, Matrix3d H, double h, double tau) {
    return Wln - tau * A * (Wl__1n - Wl_1n) / (2 * h) +
           tau * (OmegaT.inverse() * H * OmegaT) * (Wl__1n - 2 * Wln + Wl_1n) / (2 * h);
}


int main() {

    //Неизвестные
    double u;
    double e;
    //

    double L = 20.;
    double h = L / 100.;
    const int Nx = 100;
    double tau = 0.00001;
    double g = (5. / 3.);
    double Max_uc = 0;
    double Co = 0.005;
    vector<double> lol;
    double c;// speed of sound
    double t = 0.;
    double T = 0.015;
    int j = 0;
    array<Vector3d, Nx> w0;
    array<Vector3d, Nx> w1;

    // нулевой слой, начальные условия, плотность, скорость, давление
    for (int i = 0; i < Nx / 2; i++) {
        w0[i][0] = 13.;
        w0[i][1] = 0.;
        w0[i][2] = 1000000 / (g - 1);
    }
    for (int i = Nx / 2; i < Nx; i++) {
        w0[i][0] = 1.3;
        w0[i][1] = 0.;
        w0[i][2] = 100000 / (g - 1);
    }
    //
    Matrix3d OmegaT;
    Matrix3d A;
    Matrix3d H;
    while (t <= T) {
        for (int i = 0; i < Nx; i++) {
            u = w0[i][1] / w0[i][0];
            e = w0[i][2] / w0[i][0];
            c = sqrt((g - 1) * g * e);
            if (max(max(abs(u + c), u), abs(u - c)) > Max_uc) {
                Max_uc = max(max(abs(u + c), u), abs(u - c));
            }
            OmegaT << -u * c, c, g - 1,
                    -c * c, 0, g - 1,
                    u * c, -c, g - 1;
            A << 0, 1, 0,
                    -u * u, 2 * u, g - 1,
                    -g * u * e, g * e, u;
            H << abs(u + c), 0, 0,
                    0, abs(u), 0,
                    0, 0, abs(u - c);
            if ((i > 0) and (i < Nx - 1)) {
                w1[i] = KIR(w0[i - 1], w0[i], w0[i + 1], OmegaT, A, H, h, tau);
            } else if (i == 0) {
                w1[i] = KIR(w0[i], w0[i], w0[i + 1], OmegaT, A, H, h, tau);
            } else if (i == Nx - 1) {
                w1[i] = KIR(w0[i - 1], w0[i], w0[i], OmegaT, A, H, h, tau);
            }
        }
        for (int i = 0; i < Nx; i++) {
            w0[i] = w1[i];
        }
        lol.push_back(tau);
        cout << t << endl;
        t += tau;
        tau = Co * h / Max_uc;
    }
    array<double, Nx> P;
    array<double, Nx> U;
    array<double, Nx> Ro;
    array<double, Nx> E;
    for (int i = 0; i < Nx; i++) {
        P[i] = w1[i][2] * (g - 1);
        U[i] = w1[i][1] / w1[i][0];
        Ro[i] = w1[i][0];
        E[i] = w1[i][2] / w1[i][0];
    }
    ofstream out("P.txt");
    for (int x = 0; x < Nx; x++) {
        out << P[x] / 100000 << ' ';
    }
    ofstream out1("U.txt");
    for (int x = 0; x < Nx; x++) {
        out1 << U[x] << ' ';
    }
    ofstream out2("Ro.txt");
    for (int x = 0; x < Nx; x++) {
        out2 << Ro[x] << ' ';
    }
    ofstream out3("E.txt");
    for (int x = 0; x < Nx; x++) {
        out3 << E[x] << ' ';
    }
    ofstream out4("T.txt");
    for (int x = 0; x < lol.size(); x++) {
        out4 << lol[x]*1000000 << ' ';
    }
}




















