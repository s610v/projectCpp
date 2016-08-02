#include <cmath>
#include <cstdlib>
#include <iostream>
#include "IntFunction.h"
#include "gauss.hpp"
using namespace std;

#define PI 3.14159265
#define FUNC(x) ((*func)(x))
#define r0 2.81794e-13 // classical electron radius

IntFunction::IntFunction(){ // default constructor
    w1 = 1./500.;
    w2 = 1000.;
    intValues = new double*[K];
    for (int i=0; i<K; i++) {
        intValues[i] = new double[N];
    }
    epsValues = new double*[K];
    for( int i=0; i<K; i++) {
        epsValues[i] = new double[N];
    }
    w2Values = new double[K];
}
IntFunction::IntFunction(double a, double b){ // set your own w1, w2 values
    w1 = a;
    w2 = b;
    epsValues = new double*[K];
    for( int i=0; i<K; i++) {
        epsValues[i] = new double[N];
    }
    intValues = new double*[K];
    for (int i=0; i<K; i++) {
        intValues[i] = new double[N];
    }
    w2Values = new double[K];
}

double IntFunction::getW2Values(int i) {
    return w2Values[i];
}

double IntFunction::getEpsValues(int i, int j){
    return epsValues[i][j];
}

double IntFunction::getIntValues(int i, int j) {
    return intValues[i][j];
}

double IntFunction::I(double epsilon) const { // spectrum of producted particles
    double E = w2;
    double term = (E-epsilon)*epsilon;
    double sp = PI*r0*r0/(4*w1*w1*pow(w2,3)) *
     (4*E*E/term*log(4*w1*term/E) -
      8*w1*E +
      2*(2*w1*E-1)*E*E/term -
      (1-1./(w1*E))*pow(E,4)/pow(term,2));
    return sp;
}

void IntFunction::integration(){ // compute the integrated functions
    for(int m=0; m<K; m++){ // there is space left in the arrays
        w2Values[m] = w2;
        double E = w2;
        double eMin = E/2.*(1-sqrt(1-1./(w1*E))); // minimum epsilon
        double eMax = E/2.*(1+sqrt(1-1./(w1*E))); // maximum epsilon
        for(int n=0; n<N; n++){
            epsValues[m][n] = eMin + n*(eMax-eMin)/(N-1); // equal spacing of epsilon values
        }
        for(int n=0; n<N; n++){
            intValues[m][n] = qgaus<IntFunction>(*this,&IntFunction::I,eMin,epsValues[m][n]); // numerical integration to find integrated function values
            epsValues[m][n] /= eMax; // normalize epsilon
        }
        for(int n=0; n<N; n++){
            intValues[m][n] /= intValues[m][N-1]; // transform all values between 0 and 1
        }
        w2 += 100; // increases w2 to compute the next integrated function in the table
    }
}