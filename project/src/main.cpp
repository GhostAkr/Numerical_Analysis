#include "../include/PoissonEquation.h"
#include <iostream>

int main() {
    string path = "../data/outPoisson.dat";
    double t = 0.1;
    double h1 = 0.1;
    double h2 = 0.1;
    int testNum = 2;
    switchDirectionsScheme(path, t, h1, h2, testNum);
    return 0;
}