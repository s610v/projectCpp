//
//  Interpolate.cpp
//  pairproduce
//
//  Created by s610v on 7/19/16.
//  Copyright © 2016 s610v. All rights reserved.
//
#include <iostream>
#include <cmath>
#include "Interpolate.hpp"
#include "nr.h"
#include "nrutil.h"
#include "search.h"
#include "IntFunction.h"
using namespace std;

#define NR_END 1
#define FREE_ARG char*

float *Vector(long nl, long nh) // from Numerical Recipes
/* allocate a float vector with subscript range v[nl..nh] */ {
    float *v;
    v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh) // from Numerical Recipes
/* free a float vector allocated with vector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

Interpolate::Interpolate(){
    w2I = 5098.;
    cout << "interpolate object created" << endl;
}

Interpolate::Interpolate(float w2){
    w2I = w2;
    cout << "interpolate object created" << endl;
}

void Interpolate::polint(float xa[], float ya[], int n, float x, float *y, float *dy){ // From Numerical Recipes, interpolates.
    int i,m,ns=1;
    float den,dif,dift,ho,hp,w;
    float *c,*d;
    
    dif = fabs(x-xa[1]);
    c = Vector(1,n);
    d = Vector(1,n);
    for(i=1;i<=n;i++){
        if((dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif = dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1; m<n; m++) {
        for (i=1; i<=n-m; i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            den = ho-hp;
            if(den == 0.0){
                cout << "Error in routine polint" << endl;
                continue;
            }
            w=c[i+1]-d[i];
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        *y += (*dy=(2*ns < (n-m)?c[ns+1]:d[ns--]));
    }
    free_vector(d,1,n);
    free_vector(c,1,n);
}

void Interpolate::interpolation(float* fArray, Search s, IntFunction integ){ // If needed, this function gathers the necessary values from the table and interpolates
    cout << "testing" << endl;
    int isExact = s.getIsExact();
    if (isExact == 0){ // check to see if interpolation is needed
        int index = s.getIndex();
        int N = integ.getN();
        int K = integ.getK();
        float w2Val1, w2Val2, w2Val3, w2Val4;
        float iVal1[N], iVal2[N], iVal3[N], iVal4[N];
        if(index > K-4){ // special case 1
            w2Val1 = integ.getW2Values(K-4); // Isolate the correct w2 and gamma values
            w2Val2 = integ.getW2Values(K-3);
            w2Val3 = integ.getW2Values(K-2);
            w2Val4 = integ.getW2Values(K-1);
            for(int i=0; i<N; i++){
                iVal1 [i] = integ.getIntValues(K-4,i);
                iVal2 [i] = integ.getIntValues(K-3,i);
                iVal3 [i] = integ.getIntValues(K-2,i);
                iVal4 [i] = integ.getIntValues(K-1,i);
            }
        }
        else if(index < 3){ // special case 2
            w2Val1 = integ.getW2Values(0);
            w2Val2 = integ.getW2Values(1);
            w2Val3 = integ.getW2Values(2);
            w2Val4 = integ.getW2Values(3);
            for(int i=0; i<N; i++){
                iVal1 [i] = integ.getIntValues(0,i);
                iVal2 [i] = integ.getIntValues(1,i);
                iVal3 [i] = integ.getIntValues(2,i);
                iVal4 [i] = integ.getIntValues(3,i);
            }
        }
        else{ // general case
            w2Val1 = integ.getW2Values(index-1);
            w2Val2 = integ.getW2Values(index);
            w2Val3 = integ.getW2Values(index+1);
            w2Val4 = integ.getW2Values(index+2);
            for(int i=0; i<N; i++){
                iVal1 [i] = integ.getIntValues(index-1,i);
                iVal2 [i] = integ.getIntValues(index,i);
                iVal3 [i] = integ.getIntValues(index+1,i);
                iVal4 [i] = integ.getIntValues(index+2,i);
            }
        }
        cout << "test" << endl;
        // Interpolation
        for (int i=0; i<N; i++) { // loop over all epsilon
            float fw2, dfw2;
            float wVal [4] = {w2Val1, w2Val2, w2Val3, w2Val4}, iVal [4] = {iVal1[i], iVal2[i], iVal3[i], iVal4[i]};
            polint(wVal-1, iVal-1, 4, w2I, &fw2, &dfw2); // Call the polynomial interpolation class
            fArray[i] = fw2; // store interpolated values
        }
    }
}

