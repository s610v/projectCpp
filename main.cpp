//
//  main.cpp
//  pairproduce
//
//  Created by s610v on 7/18/16.
//  Copyright © 2016 s610v. All rights reserved.
//

#include <iostream>
#include <random>
#include <ctime>
#include <fstream>
#include "search.h"
#include "IntFunction.h"
#include "Interpolate.hpp"
using namespace std;

// Cutpoint method
void CM(int index, float Z[], float F[], IntFunction t){ // Cutpoint method
    int N = t.getN();
    float eValues[N];
    for (int i=0; i<N; i++) {
        eValues[i] = t.getIntValues(index,0,i);
    }
    int m = N-2;
    int p[m+1];
    for (int j=0; j<m; j++) {
        for (int i=1; i<N; i++) {
            if(F[i]>(float(j)/m)) {
                p[j] = eValues[i];
                break;
            }
        }
    }
    p[N-1] = 1.0;
    for (int i=0; i<N; i++) {
        double U = (double) rand()/RAND_MAX;
        int L = int(ceil(m*U));
        int a = int(p[L]);
        while (U > F[a]) {
            a++;
        }
        Z[i] = eValues[a];
    }
}

int main(int argc, const char * argv[]) {
    srand(time(NULL));
    ofstream f;
    f.open("histdata.txt"); // Stores values for histogram
    float w2I = 8439.; // user input w2 value
    
    IntFunction integ; // integration object
    Search s(w2I); // searching object
    integ.integration(); // compute the table of integrated functions
    float w2Values[integ.getK()]; // for keeping the stored w2 values handy
    for (int i=0; i<integ.getN(); i++) {
        w2Values[i] = integ.getW2Values(i);
    }
    s.binarySearch(w2Values, integ.getK()); // search the table
    float fArray[integ.getN()];
    if(s.getIsExact() == 0){ // we must interpolate
        Interpolate intpl; // interpolation object
        intpl.interpolation(fArray, w2I, s, integ);
    }
    else if(s.getIsExact() == 1){ // we do not need interpolation
        for (int i=0; i<integ.getN(); i++) {
            fArray[i] = integ.getIntValues(s.getIndex(), 1, i);
        }
    }
    float Z[integ.getN()];
    CM(s.getIndex(),Z,fArray,integ); // go through the cutpoint method
    for (int i=0; i<integ.getN(); i++) {
        f << Z[i] << endl; // store the results in a file
    }
    f.close();
    
    ofstream g; // Parameters needed for graphics
    g.open("parameter.txt");
    g << w2I << endl;
    g << integ.getW1() << endl;
    g << integ.getN() << endl;
    g.close();
    
    cout << "Success, Project has reached its conclusion\n";
    return 0;
}