#include "../include/PoissonEquation.h"
#include <iostream>

int main() {
    string path = "../data/outPoisson.dat";
    double t = 0.05;
    double h1 = 0.05;
    double h2 = 0.05;
    cout << "h1 = " << h1 << endl;
    cout << "h2 = " << h2 << endl;
    int testNum = 2;
    switchDirectionsScheme(path, t, h1, h2, testNum);
    return 0;
}