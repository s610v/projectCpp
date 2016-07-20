#include <cmath>
#include <cstdlib>
#include "nr.h"
#include "IntFunction.h"
using namespace std;

#define PI 3.14159265
#define EPS 1.0e-6
#define JMAX 20
#define FUNC(x) ((*func)(x))

IntFunction::IntFunction(){
    w1 = 1./500.;
    w2 = 1000.;
}

IntFunction::IntFunction(float a, float b){
    w1 = a;
    w2 = b;
}

float IntFunction::I(float epsilon){
    float E = w2;
    float i = PI*r0*r0/(4*w1*w1*pow(w2,3))*(4*E*E/((E-epsilon)*epsilon)*log(4*w1*(E-epsilon)*epsilon/E)-8*w1*E+2*(2*w1*E-1)*E*E/((E-epsilon)*epsilon)-(1-1./(w1*E))*pow(E,4)/((E-epsilon)*(E-epsilon)*epsilon*epsilon));
    return i;
}

float IntFunction::qgaus(float (*func) (float), float a, float b){
    int j;
    float xr, xm, dx, s;
    static float x[] = {0.0,0.1488743389,0.4333953941,0.6794095682,0.8650633666,0.9739065285};
    static float w[] = {0.0,0.2955242247,0.2692667193,0.2190863625,0.1494513491,0.0666713443};
    
    xm=0.5*(b+a);
    xr=0.5*(b-a);
    s=0;
    for (j=1; j<=5; j++) {
        dx=xr*x[j];
        s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
    }
    return s *= xr;
}

void IntFunction::integration(){ // compute the integrated functions
    int m = 0;
    double epsilon[N];
    while(m<K){
        w2Values[m] = w2;
        float E = w2;
        float eMin = E/2.*(1-sqrt(1-1./(w1*E)));
        float eMax = E/2.*(1+sqrt(1-1./(w1*E)));
        for(int n=0; n<N; n++){
            epsilon[n] = eMin + n*(eMax-eMin)/(N-1); // equal spacing for N values
        }
        float iArray [N];
        for(int n=0; n<N; n++){
            iArray[n] = qgaus(&I,eMin,epsilon[n]);
            epsilon[n] /= eMax;
            intValues[m][0][n] = epsilon[n];
        }
        for(int n=0; n<N; n++){
            iArray[n] /= iArray[N-1]; // transform all values between 0 and 1
            intValues[m][1][n] = iArray[n];
        }
        w2 += 100;
        m++;
    }
}
