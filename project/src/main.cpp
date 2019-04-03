#include "../include/OrdinaryDE.h"
#include "../include/HeatEquation.h"
#include <iostream>

int main() {
    string path = "../data/out.dat";
    mixedScheme(0.5, 0.02, 0.0002, path);
    return 0;
}