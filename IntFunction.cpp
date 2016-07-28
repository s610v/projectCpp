#include <cmath>
#include <cstdlib>
#include <iostream>
#include "IntFunction.h"
#include "gauss.hpp"
using namespace std;
#define PI 3.14159265
#define FUNC(x) ((*func)(x))

IntFunction::IntFunction(){
    w1 = 1./500.;
    w2 = 1000.;
}
IntFunction::IntFunction(float a, float b){
    w1 = a;
    w2 = b;
}

const float IntFunction::r0 = 2.81794e-13;

float IntFunction::I(float epsilon) const {
    float E = w2;
    float sp = PI*r0*r0/(4*w1*w1*pow(w2,3))*(4*E*E/((E-epsilon)*epsilon)*log(4*w1*(E-epsilon)*epsilon/E)-8*w1*E+2*(2*w1*E-1)*E*E/((E-epsilon)*epsilon)-(1-1./(w1*E))*pow(E,4)/((E-epsilon)*(E-epsilon)*epsilon*epsilon));
    return sp;
}

void IntFunction::integration(){ // compute the integrated functions
    int m = 0;
    while(m<K){
        w2Values[m] = w2;
        float E = w2;
        float eMin = E/2.*(1-sqrt(1-1./(w1*E)));
        float eMax = E/2.*(1+sqrt(1-1./(w1*E)));
        for(int n=0; n<N; n++){
            epsValues[n] = eMin + n*(eMax-eMin)/(N-1); // equal spacing for N values
        }
        for(int n=0; n<N; n++){
            intValues[m][n] = qgaus<IntFunction>(*this,&IntFunction::I,eMin,epsValues[n]);
            epsValues[n] /= eMax;
        }
        for(int n=0; n<N; n++){
            intValues[m][n] /= intValues[m][N-1]; // transform all values between 0 and 1
        }
        w2 += 100;
        m++;
    }
}