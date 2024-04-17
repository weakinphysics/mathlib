#include "./Ultimatrix.h"
#include <iostream>

using namespace std;

int main(){
    double p[9] = {1,2,3,4,5,6,7,8,9};
    double q[9] = {9,8,7,6,5,4,3,2,1};
    Ultimatrix <double> e(3, 3, p);
    e.setElement(0, 0, 1.00);
    e.setElement(0, 1, 1.00);
    e.setElement(0, 2, 1.00);
    e.setElement(1, 0, 1.00);
    e.setElement(2, 0, 1.00);
    e.setElement(1, 1, 2.00);
    e.setElement(1, 2, 2.00);
    e.setElement(2, 1, 2.00);
    e.setElement(2, 2, 3.00);    
    // Ultimatrix<double> inv = e.invert();
    vector<Ultimatrix<double>> qr = e.qrf();
    qr[0].printMatrix();
    qr[1].printMatrix();
    // inv.printMatrix();
    return 0;
}