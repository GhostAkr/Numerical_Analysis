#include "../include/WaveEquation.h"
#include <iostream>

int main() {
    string path = "../data/waveOut.dat";
    double t = 0.000625;
    double h = 0.005;
    crossScheme(t, h, path, 4);
    return 0;
}