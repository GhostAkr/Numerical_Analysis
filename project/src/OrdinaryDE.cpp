//
// Created by ighos on 08.02.2019.
//

#include "../include/OrdinaryDE.h"

vector<double> func0(vector<double> _point) {
    // Constants
    double k = 20.0;
    double m = 0.3;
    if (_point.size() != 2) {  // Exception
        cout << "Incorrect point" << endl;
        return {};
    }
    vector<double> result(2);
    result[0] = _point[1];
    result[1] = -k / m * _point[0];
    return result;
}

void vectorPrint(vector<double> _sourceVector) {
    for (int i = 0; i < _sourceVector.size(); ++i) {
        cout << _sourceVector[i] << endl;
    }
    cout << endl;
}

void EulerExplicit(vector<double> _func(vector<double>), vector<double> _startPoint) {
    std::ofstream fOut("../data/solution.dat");
    if (!fOut) {  // Exception
        cout << "Error while opening file" << endl;
        return;
    }
    double step = 0.01;
    double maxMesh = 500;
    //fOut << maxMesh << endl;
    size_t nOfVars = _startPoint.size();
    //cout << nOfVars << endl;
    for (int i = 1; i < maxMesh; ++i) {
        vector<double> point(nOfVars);
        vector<double> f = _func(_startPoint);
        for (int j = 0; j < nOfVars; ++j) {
            point[j] = _startPoint[j] + step * f[j];
            fOut << point[j] << " ";
        }
        fOut << endl;
        _startPoint = point;
    }
}
