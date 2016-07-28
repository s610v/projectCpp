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
#include <string>
#include <cstdlib>
#include "search.h"
#include "IntFunction.h"
#include "Interpolate.hpp"
using namespace std;

// Cutpoint method
void CM(float F[], IntFunction integ){ // Cutpoint method
    int N = integ.getN();
    int m = N-2;
    int p[m+1];
    for (int j=0; j<m; j++) {
        for (int i=1; i<N; i++) {
            if(F[i]>(float(j)/m)) {
                p[j] = integ.getEpsValues(i);
                break;
            }
        }
    }
    p[m] = 1.0;
    ofstream f;
    f.open("histdata.txt"); // Stores values for histogram
    for (int i=0; i<N; i++) {
        float U = (float) rand()/RAND_MAX;
        int L = int(ceil(m*U));
        int a = int(p[L]);
        while (U > F[a]) {
            a++;
        }
        f << integ.getEpsValues(a) << endl; // store the results in a file
    }
    f.close();
}

int main(int argc, const char * argv[]) {
    float w1 = 1./500.; // energy of soft photons
    float w2 = 1000.; // initial energy value of gamma ray
    float w2I;
    cout << "Enter a float number: ";
    cin >> w2I; // users enter a float for input w2
    IntFunction integ(w1, w2); // integration object
    Search search(w2I); // searching object
    srand(1);
    integ.integration(); // compute the table of integrated functions
    //float* w2Values = new float[integ.getK()];
    //float* w2Values = (float*) malloc(integ.getK();
    float w2Values[integ.getK()]; // for storing the w2 values
    for (int i=0; i<integ.getK(); i++) {
        w2Values[i] = integ.getW2Values(i);
    }
    search.binarySearch(w2Values, integ.getK()); // search the table
    //delete[] w2Values;
    //free(w2Values);
    float fArray[integ.getN()];
    //float* fArray = new float[integ.getN()];
    if(search.getIsExact() == 0){ // we must interpolate
        Interpolate intpl(w2I); // interpolation object
        intpl.interpolation(fArray, search, integ);
        CM(fArray,integ); // go through the cutpoint method
    }
    else if(search.getIsExact() == 1){ // we do not need interpolation
        for (int i=0; i<integ.getN(); i++) {
            fArray[i] = integ.getIntValues(search.getIndex(), i);
        }
        CM(fArray,integ); // go through the cutpoint method
    }
    //delete[] fArray;
    ofstream g; // Parameters needed for graphics
    g.open("parameter.txt");
    g << w2I << endl;
    g << integ.getW1() << endl;
    g << integ.getN() << endl;
    g.close();
    
    cout << "Success, Project has reached its conclusion.\n";
    return 0;
}