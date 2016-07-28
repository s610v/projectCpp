// Header file for IntFunction.cpp

#ifndef IntFunction_h
#define IntFunction_h

#include <stdio.h>

class IntFunction{
    float w1; // fixed energy of soft photons
    float w2; // variable energy of gamma rays
    static const float r0; // Classical electron radius
    static const int N = 10000; // number of points for arrays
    static const int K = 70; // number of w2 values used
    float intValues[K][N];
    float epsValues [N];
    float w2Values [K];
public:
    IntFunction(); // Default constructor
    IntFunction(float, float); // set your own values
    float I(float) const; // Spectrum function
    void integration(); // calls the Gauss quadrature method, also sets up integration
    float getIntValues(int i, int k) {return intValues[i][k];} // returns the array with epsilon, integrated function values
    float getEpsValues(int i) {return epsValues[i];}
    float getW2Values(int i) {return w2Values[i];} // returns the array of w2 values
    int getK() {return K;} // returns the number of w2 values used
    int getN() {return N;} // returns the number of epsilon values used
    float getW1() {return w1;} // returns the fixed value of w1
    float getW2() {return w2;}
};
#endif /* IntFunction_h */