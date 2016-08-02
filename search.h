// Header file for search.cpp


#ifndef search_h
#define search_h

#include <stdio.h>
class Search{
    int index; // to store the index closest to the input w2 value
    int isExact; // to store whether or not the input value is on the computed table
    double w2In; // user input w2 value
public:
    Search();
    Search(double); // sets w2 input value
    int getIndex() {return index;} // returns index of closest value in table
    int getIsExact() {return isExact;} // returns isExact
    double getW2In() {return w2In;} // returns your input value
    void binarySearch(double*, int);
};
#endif /* search_h */