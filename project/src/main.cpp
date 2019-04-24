#include "../include/WaveEquation.h"
#include <iostream>

int main() {
    string path = "../data/waveOut.dat";
    double t = 0.025;
    double h = 0.05;
    crossScheme(t, h, path, 2);
    return 0;
}