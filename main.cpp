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
#include "Interpolate.h"
using namespace std;

// Cutpoint method
void CM(double F[], IntFunction integ, int index){
    int N = integ.getN();
    int m = N-2; // set the number of pointers
    int* p = new int[m+1]; // stores the pointers
    for (int j=0; j<m; j++) { // compute the pointers
        for (int i=1; i<N; i++) {
            if(F[i]>(double(j)/m)) {
                p[j] = integ.getEpsValues(index,i);
                break;
            }
        }
    }
    p[m] = 1.0; // last pointer is the upper bound of normalized epsilon, which is always 1.0
    ofstream f;
    f.open("histdata.txt"); // Stores values for histogram
    for (int i=0; i<N; i++) {
        double U = (double) rand()/RAND_MAX; // get random number
        int L = int(ceil(m*U));
        int a = int(p[L]);
        while (U > F[a]) {
            a++;
        }
        f << integ.getEpsValues(index,a) << endl; // store the results in a file
    }
    delete[] p;
    f.close();
}

int main() {
    IntFunction integ(1./500., 1000.); // integration object, w1 and w2 are the arguments
    double w2I; // input w2 Value
    cout << "Enter a number ";
    cin >> w2I;
    if(w2I < integ.getW2()){ // catch some possible out of bounds early
        cout << "w2 input value is below the minimum of " << integ.getW2() << " (or is of invalid type), ending program. Try again." << endl;
        return 0;
    }
    Search search(w2I); // searching object
    integ.integration(); // compute the table of integrated functions
    cout << "Integration Complete" << endl;
    
    double* w2Values = new double[integ.getK()];
    for (int i=0; i<integ.getK(); i++) {
        w2Values[i] = integ.getW2Values(i); // store w2 values locally
    }
    search.binarySearch(w2Values, integ.getK()); // search the table
    delete[] w2Values;
    
    if(search.getIsExact() == -1){
        cout << "Program has ended due to out of bounds value." << endl;
        return 0;
    }
    
    double* fArray = new double[integ.getN()];
    if(search.getIsExact() == 0){ // we must interpolate
        Interpolate intpl;
        intpl.interpolation(fArray, w2I, search, integ); // calls interpolation function
        integ.delW2Values(); // deletes the array of w2Values to free up stack space
    }
    else if(search.getIsExact() == 1){ // we do not need interpolation
        integ.delW2Values();
        for (int i=0; i<integ.getN(); i++) {
            fArray[i] = integ.getIntValues(search.getIndex(), i); // take the values from the table
        }
    }
    
    srand(time(NULL)); // set the seed of the random number
    CM(fArray,integ, search.getIndex()); // go through the cutpoint method
    delete[] fArray;
    ofstream g; // Parameters needed for graphics
    g.open("parameter.txt");
    g << w2I << endl;
    g << integ.getW1() << endl;
    g << integ.getN() << endl;
    g.close();
    
    cout << "Success, project has reached its conclusion.\n";
    return 0;
}