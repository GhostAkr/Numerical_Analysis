#include "../include/OrdinaryDE.h"
#include "../include/HeatEquation.h"
#include <iostream>

int main() {
    string path = "../data/out.dat";
    mixedScheme(1, 0.01, 0.01, path);
    //quasilinearEquationScheme(0.2, 2e-4, path);
    //exactQuasiChart(0.2, 2e-4, path);
    return 0;
}