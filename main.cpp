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
    double aadf = 10e-6

}




















//double
//KIR1(double Wl_1n1, double Wln1, double Wl__1n1, double Wl_1n2, double Wln2, double Wl__1n2, double Wl_1n3, double Wln3,
//     double Wl__1n3, double tau, double h, double A1, double A2, double A3, double H1,
//     double H2, double H3) {
//    return Wln1 - tau / (2 * h) * (A1 * (Wl__1n1 - Wl_1n1) + A2 * (Wl__1n2 - Wl_1n2) + A3 * (Wl__1n3 - Wl_1n3)) +
//           tau / (2 * h) * (H1 * (Wl__1n1 - 2 * Wln1 + Wl_1n1) + H2 * (Wl__1n2 - 2 * Wln2 + Wl_1n2) +
//                            H3 * (Wl__1n3 - 2 * Wln3 + Wl_1n3));
//}
//
//double
//KIR2(double Wl_1n1, double Wln1, double Wl__1n1, double Wl_1n2, double Wln2, double Wl__1n2, double Wl_1n3, double Wln3,
//     double Wl__1n3, double tau, double h, double A1, double A2, double A3, double H1,
//     double H2, double H3) {
//    return Wln2 - tau / (2 * h) * (A1 * (Wl__1n1 - Wl_1n1) + A2 * (Wl__1n2 - Wl_1n2) + A3 * (Wl__1n3 - Wl_1n3)) +
//           tau / (2 * h) * (H1 * (Wl__1n1 - 2 * Wln1 + Wl_1n1) + H2 * (Wl__1n2 - 2 * Wln2 + Wl_1n2) +
//                            H3 * (Wl__1n3 - 2 * Wln3 + Wl_1n3));
//}
//
//double
//KIR3(double Wl_1n1, double Wln1, double Wl__1n1, double Wl_1n2, double Wln2, double Wl__1n2, double Wl_1n3, double Wln3,
//     double Wl__1n3, double tau, double h, double A1, double A2, double A3, double H1,
//     double H2, double H3) {
//    return Wln3 - tau / (2 * h) * (A1 * (Wl__1n1 - Wl_1n1) + A2 * (Wl__1n2 - Wl_1n2) + A3 * (Wl__1n3 - Wl_1n3)) +
//           tau / (2 * h) * (H1 * (Wl__1n1 - 2 * Wln1 + Wl_1n1) + H2 * (Wl__1n2 - 2 * Wln2 + Wl_1n2) +
//                            H3 * (Wl__1n3 - 2 * Wln3 + Wl_1n3));
//}
//
//
//int main() {
//    double L = 20.;
//    double h = L / 100.;
//    int Nx = 100;
//    double tau = 0.000001;
//    double g = (5. / 3.);
//    double Co = 0.005;
//    vector<vector<double>> w1;
//    w1.assign(2, vector<double>(Nx));// ro
//    vector<vector<double>> w2;
//    w2.assign(2, vector<double>(Nx));// ro*u
//    vector<vector<double>> w3;
//    w3.assign(2, vector<double>(Nx));// ro*e
//    vector<vector<double>> A;
//    vector<vector<double>> H;
//    A.assign(3, vector<double>(3));
//    H.assign(3, vector<double>(3));
//    //Начальные условия >>
//    for (int i = 0; i < Nx / 2; i++) {
//        w1[0][i] = 13.;
//        w2[0][i] = 0.;
//        w3[0][i] = 1000000 / (g - 1);
//    }
//    for (int i = Nx / 2; i < Nx; i++) {
//        w1[0][i] = 1.3;
//        w2[0][i] = 0.;
//        w3[0][i] = 100000 / (g - 1);
//    }
//    A[0][0] = 0.;
//    A[0][1] = 1.;
//    A[0][2] = 0.;
//    A[1][2] = g - 1.;
//    double e;
//    double c;
//    double u;
//    double t = 0;
//    double T = 0.05;
//
//    //<<
//    t += tau;
//    while (t <= T) {
//        for (int i = 0; i < Nx; i++) {
//            A[1][0] = -(w2[0][i] / w1[0][i]) * (w2[0][i] / w1[0][i]);
//            A[1][1] = w2[0][i] / w1[0][i] * 2.;
//            A[2][0] = -g * (w2[0][i] / w1[0][i]) * (w3[0][i] / w1[0][i]);
//            A[2][1] = g * (w3[0][i] / w1[0][i]);
//            A[2][2] = w2[0][i] / w1[0][i];
//            e = w3[0][i] / w1[0][i];
//            c = sqrt(g * (g - 1) * e);
//            u = A[2][2];
//            double a_1 = abs(u + c);
//            double a1 = abs(u - c);
//            double a = abs(u);
//            //lutiy pizdec, multiplied matrix before last part of equation, looks like shit and i should remake this
//            H[0][0] = -u * a_1 / (2. * c) + a + a1 * u / (2. * c);
//            H[0][1] = (a_1 - a1) / (2. * c);
//            H[0][2] = (g - 1.) * (a_1 / (2. * c * c) - a / (c * c) + a1 / (2. * c * c));
//            H[1][0] = -a_1 * (u + c) * u / (2. * c) + u * a + a1 * (u - c) * u / (2. * c);
//            H[1][1] = a_1 * (u + c) / (2. * c) - a1 * (u - c) / (2. * c);
//            H[1][2] = (g-1)*(c*(a_1-a1)+u*(a_1+a1) - 2*u*a)/(2*c*c);
//            H[2][0] = (a_1 - a1) * u * c / (2. * (g - 1.));
//            H[2][1] = (a_1 - a1) * c / (2. * (g - 1.));
//            H[2][2] = (a_1 + a1) / 2.;
//            if ((i != 0) and i != (Nx - 1)) {
//                w1[1][i] = KIR1(w1[0][i - 1], w1[0][i], w1[0][i + 1], w2[0][i - 1], w2[0][i], w2[0][i + 1],
//                                w3[0][i - 1],
//                                w3[0][i], w3[0][i + 1], tau, h, A[0][0], A[0][1], A[0][2], H[0][0],
//                                H[0][1], H[0][2]);
//                w2[1][i] = KIR2(w1[0][i - 1], w1[0][i], w1[0][i + 1], w2[0][i - 1], w2[0][i], w2[0][i + 1],
//                                w3[0][i - 1],
//                                w3[0][i], w3[0][i + 1], tau, h, A[1][0], A[1][1], A[1][2], H[1][0],
//                                H[1][1], H[1][2]);
//                w3[1][i] = KIR3(w1[0][i - 1], w1[0][i], w1[0][i + 1], w2[0][i - 1], w2[0][i], w2[0][i + 1],
//                                w3[0][i - 1],
//                                w3[0][i], w3[0][i + 1], tau, h, A[2][0], A[2][1], A[2][2], H[2][0],
//                                H[2][1], H[2][2]);
//            }
//            else if(i == 0){
//                w1[1][i] = KIR1(w1[0][i], w1[0][i], w1[0][i + 1], w2[0][i], w2[0][i], w2[0][i + 1],
//                                w3[0][i],
//                                w3[0][i], w3[0][i + 1], tau, h, A[0][0], A[0][1], A[0][2], H[0][0],
//                                H[0][1], H[0][2]);
//                w2[1][i] = KIR2(w1[0][i], w1[0][i], w1[0][i + 1], w2[0][i], w2[0][i], w2[0][i + 1],
//                                w3[0][i],
//                                w3[0][i], w3[0][i + 1], tau, h, A[1][0], A[1][1], A[1][2], H[1][0],
//                                H[1][1], H[1][2]);
//                w3[1][i] = KIR3(w1[0][i], w1[0][i], w1[0][i + 1], w2[0][i], w2[0][i], w2[0][i + 1],
//                                w3[0][i],
//                                w3[0][i], w3[0][i + 1], tau, h, A[2][0], A[2][1], A[2][2], H[2][0],
//                                H[2][1], H[2][2]);
//            }
//            else if(i == (Nx-1)){
//                w1[1][i] = KIR1(w1[0][i - 1], w1[0][i], w1[0][i], w2[0][i - 1], w2[0][i], w2[0][i],
//                                w3[0][i - 1],
//                                w3[0][i], w3[0][i], tau, h, A[0][0], A[0][1], A[0][2], H[0][0],
//                                H[0][1], H[0][2]);
//                w2[1][i] = KIR2(w1[0][i - 1], w1[0][i], w1[0][i], w2[0][i - 1], w2[0][i], w2[0][i],
//                                w3[0][i - 1],
//                                w3[0][i], w3[0][i], tau, h, A[1][0], A[1][1], A[1][2], H[1][0],
//                                H[1][1], H[1][2]);
//                w3[1][i] = KIR3(w1[0][i - 1], w1[0][i], w1[0][i], w2[0][i - 1], w2[0][i], w2[0][i],
//                                w3[0][i - 1],
//                                w3[0][i], w3[0][i], tau, h, A[2][0], A[2][1], A[2][2], H[2][0],
//                                H[2][1], H[2][2]);
//            }
//
//        }
//        for (int i = 0; i < Nx; i++) {
//            w1[0][i] = w1[1][i];
//            w2[0][i] = w2[1][i];
//            w3[0][i] = w3[1][i];
//        }
//        tau = Co * h / max(max(abs(u + c), abs(u - c)), abs(u));
//        t += tau;
//        cout << t << endl;
//    }
//    //непосредственно вывод
//    vector<double> P;
//    P.resize(Nx);
//    vector<double> Ro;
//    Ro.resize(Nx);
//    vector<double> U;
//    U.resize(Nx);
//    vector<double> E;
//    E.resize(Nx);
//    for (int i = 0; i < Nx; i++) {
//        P[i] = (g - 1) * w3[0][i];
//        Ro[i] = w1[1][i];
//        U[i] = w2[1][i] / w1[1][i];
//        E[i] = w3[1][i] / w1[1][i];
//    }
//    ofstream out("U.txt");
//    for (int x = 0; x < Nx; x++) {
//        out << U[x] << ' ';
//    }
//}



