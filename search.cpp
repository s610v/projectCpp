#include <cstdlib>
#include <iostream>
#include <string>
#include "search.h"
using namespace std;

Search::Search(){
    index = isExact = 0;
    w2In = 5098.; // user input w2 value
}

Search::Search(double w2){ // allows the user to set their own w2 value for input
    index = isExact = 0; // Both are preset to zero
    w2In = w2;
}

void Search::binarySearch(double* w2Values, int K){ // searches the table to find closest values to input value
    int a = 0;
    int b = K - 1;
    while((b-a)>1) {
        if((w2In < w2Values[a]) || w2In > w2Values[b]){ // w2In was not in bounds of the table's values
            cout << "Invalid w2 input. Try again." << endl;
            isExact = -1;
            break;
        }
        else if(w2In == w2Values[a]){ // w2In was the lower bound
            isExact = 1;
            index = a;
            break;
        }
        else if(w2In == w2Values[b]){ // w2In was the upper bound
            isExact = 1;
            index = b;
            break;
        }
        
        int mid = int(double(b+a)/2.); // calculate the midpoint of the interval
        
        if (w2In > w2Values[mid]) { // update the interval
            a = mid;
            index = b;
        }
        else if(w2In < w2Values[mid]) { // update the interval
            b = mid;
            index = a;
        }
        else { // input was the midpoint's value
            index = mid;
            isExact = 1;
            break;
        }
    }
    cout << "Table search complete. " << endl;
}