#include "../include/OrdinaryDE.h"
#include "../include/HeatEquation.h"
#include <iostream>

int main() {
    string path = "../data/out.dat";
    mixedScheme(1, 0.1, 0.001, path);
    return 0;
}