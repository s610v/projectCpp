// Header file for IntFunction.cpp

#ifndef IntFunction_h
#define IntFunction_h

#include <stdio.h>
#include <cmath>

class IntFunction{
    float w1; // fixed energy of soft photons
    float w2; // variable energy of gamma rays
    constexpr static const float r0 = 2.81794E-13; // Classical electron radius
    static const int N = 600; // number of points for arrays
    static const int K = 76; // number of w2 values used
    float intValues [K][2][N];
    float w2Values [K];
public:
    IntFunction(); // Default constructor
    IntFunction(float, float); // set your own values
    float I(float); // Spectrum function
    void integration(); // calls the Gauss quadrature method, also sets up integration
    float getIntValues(int i, int j, int k) {return intValues[i][j][k];} // returns the array with epsilon, integrated function values
    float getW2Values(int i) {return w2Values[i];} // returns the array of w2 values
    int getK() {return K;} // returns the number of w2 values used
    int getN() {return N;} // returns the number of epsilon values used
    float getW1() {return w1;} // returns the fixed value of w1
    float getR0() {return r0;} // returns r0
    float qgaus(float (*func)(float), float, float);
};
#endif /* IntFunction_h */