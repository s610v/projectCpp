// Header file for IntFunction.cpp

#ifndef IntFunction_h
#define IntFunction_h

#include <stdio.h>

class IntFunction{
    double w1; // fixed energy of soft photons
    double w2; // variable energy of gamma rays
    static const int N = 10000; // number of points for arrays
    static const int K = 76; // number of w2 values used
    double** intValues; // stores integrated function values
    double** epsValues; // stores normalized epsilon values
    double* w2Values; // stores the w2 values from the table
public:
    IntFunction(); // Default constructor
    IntFunction(double, double); // set your own values
    double I(double) const; // Spectrum function
    void integration(); // calls the Gauss quadrature method, also sets up integration
    double getIntValues(int, int); // returns a value from an array with integrated function values
    double getEpsValues(int, int); // returns a value of epsilon from the stored array
    double getW2Values(int); // returns a w2 value from the stored array
    void delW2Values() {delete[] w2Values;} // deletes the array of w2 values when it is not needed anymore
    int getK() const {return K;} // returns the number of w2 values used
    int getN() const {return N;} // returns the number of epsilon values used
    double getW1() {return w1;} // returns the fixed value of w1
    double getW2() {return w2;} // returns w2
};
#endif /* IntFunction_h */