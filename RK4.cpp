#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <cstdlib>

using namespace std;

#define PI 3.1415926535897932
#define n 2

double t0 = 0;
double tmax = 20;
double dt = 0.01;

double om = 2*PI;
double om_0 = 1.5*om;
double b = om_0 * 0.25;
double g = 1.073;

double ics[n] = {0, 1};

// Defines you differential equation
double f(double t, double x[n], int counter){
    if (counter == 0) {
        return x[1];
    } else if (counter == 1) {
        return -2*b*x[1] - (om_0*om_0)*sin(x[0]) + g*(om_0*om_0)*cos(om*t);
    };
};

int main() {
    double t;
    double x[n], k1[n], k2[n], k3[n], k4[n];
    double tempk1[n], tempk2[n], tempk3[n];

    t = t0;
    copy(ics +0, ics +n, x);

    ofstream myfile("Rk4_sol.dat") ;
    myfile.precision (17);

    while (t <= tmax){
        myfile << t << " ";
        for (int i = 0; i < n; i++) {
            myfile << x[i] << " ";
        };
        myfile << endl;

        for (int i = 0; i < n; i++) {
            k1[i] = dt * f(t, x, i);
            tempk1[i] = x[i] + 0.5*k1[i];
        };
        for (int i = 0; i < n; i++) {
            k2[i] = dt * f(t + 0.5*dt, tempk1, i);
            tempk2[i] = x[i] + 0.5*k2[i];
        };
        for (int i = 0; i < n; i++) {
            k3[i] = dt * f(t + 0.5*dt, tempk2, i);
            tempk3[i] = x[i] + k3[i];
        };
        for (int i = 0; i < n; i++) {
            k4[i] = dt * f(t + dt, tempk3, i);
        };

        for (int i = 0; i < n; i++) {
            x[i] += (1.0/6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        };
        t += dt;
    };

    return 0;
}