//
//  Interpolate.cpp
//  pairproduce
//
//  Created by s610v on 7/19/16.
//  Copyright © 2016 s610v. All rights reserved.
//
#include <iostream>
#include "Interpolate.hpp"
#include "nr.h"
#include "nrutil.h"
#include "search.h"
#include "IntFunction.h"
using namespace std;

void Interpolate::polint(float xa[], float ya[], int n, float x, float *y, float *dy){
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

void Interpolate::interpolation(float* fArray, float w2, Search s, IntFunction t){
    int isExact = s.getIsExact();
    int index = s.getIndex();
    int N = t.getN();
    int K = t.getK();
    if (isExact == 0){
        float w2Val1, w2Val2, w2Val3, w2Val4;
        float iVal1[t.getN()], iVal2[t.getN()], iVal3[t.getN()], iVal4[t.getN()];
        if(index > K-4){
            w2Val1 = t.getW2Values(K-4); // Isolate the correct w2 and gamma values
            w2Val2 = t.getW2Values(K-3);
            w2Val3 = t.getW2Values(K-2);
            w2Val4 = t.getW2Values(K-1);
            for(int i=0; i<t.getN(); i++){
                iVal1 [i] = t.getIntValues(K-4,1,i);
                iVal2 [i] = t.getIntValues(K-3,1,i);
                iVal3 [i] = t.getIntValues(K-2,1,i);
                iVal4 [i] = t.getIntValues(K-1,1,i);
            }
        }
        else if(index < 3){
            w2Val1 = t.getW2Values(0);
            w2Val2 = t.getW2Values(1);
            w2Val3 = t.getW2Values(2);
            w2Val4 = t.getW2Values(3);
            for(int i=0; i<t.getN(); i++){
                iVal1 [i] = t.getIntValues(0,1,i);
                iVal2 [i] = t.getIntValues(1,1,i);
                iVal3 [i] = t.getIntValues(2,1,i);
                iVal4 [i] = t.getIntValues(3,1,i);
            }
        }
        else{
            w2Val1 = t.getW2Values(index-1);
            w2Val2 = t.getW2Values(index);
            w2Val3 = t.getW2Values(index+1);
            w2Val4 = t.getW2Values(index+2);
            for(int i=0; i<t.getN(); i++){
                iVal1 [i] = t.getIntValues(index-1,1,i);
                iVal2 [i] = t.getIntValues(index,1,i);
                iVal3 [i] = t.getIntValues(index+1,1,i);
                iVal4 [i] = t.getIntValues(index+2,1,i);
            }
        }
        
        // Interpolation
        for (int i=0; i<N; i++) { // loop over all epsilon
            float fw2, dfw2;
            float wVal [4] = {w2Val1, w2Val2, w2Val3, w2Val4}, iVal [4] = {iVal1[i], iVal2[i], iVal3[i], iVal4[i]};
            polint(wVal-1, iVal-1, 4, w2, &fw2, &dfw2); // Call the polynomial interpolation class
            cout << w2 << " " << fw2 << " " << dfw2 << endl;
            fArray[i] = fw2; // store interpolated values
        }
    }
}

