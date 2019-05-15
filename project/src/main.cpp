#include "../include/WaveEquation.h"
#include <iostream>

int main() {
    string path = "../data/waveOut.dat";
    double t = 0.005;
    double h = 0.05;
//    crossScheme(t, h, path, 1);
    crossSchemeDerivative(t, h, path, 1);
    return 0;
}