#include <cstdlib>
#include <iostream>
#include <string>
#include "search.h"
using namespace std;

Search::Search(){
    index = isExact = 0;
    w2In = 8439.; // user input w2 value
}

Search::Search(float w2){ // allows the user to set their own w2 value for input
    index = isExact = 0; // Both are preset to zero always
    w2In = w2;
}

void Search::binarySearch(float* w2Values, int K){ // binary searches the table to find the closest values to input
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
        int mid = int(float(b+a)/2.); // calculate the midpoint of the interval
        if (w2In > w2Values[mid]) { // update the interval
            a = mid;
            index = b;
        }
        else if(w2In < w2Values[mid]) { // update the interval
            b = mid;
            index = a;
        }
        else{ // input was the midpoint's value
            index = mid;
            isExact = 1;
            break;
        }
    }
    cout << "Your w2 value was near the w2 value in the table with index " << index << endl;
    cout << isExact << endl; // Shows the user if their w2 value was from the table or not, 1 means that it was
}