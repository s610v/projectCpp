// Header file for search.cpp


#ifndef search_h
#define search_h

#include <stdio.h>
class Search{
    int index; // to store the index closest to the input w2 value
    int isExact; // to store whether or not the input value is on the computed table
    float w2In; // user input w2 value
public:
    Search();
    Search(float);
    int getIndex() {return index;}
    int getIsExact() {return isExact;}
    double getW2In() {return w2In;}
    void binarySearch(float*, int);
};
#endif /* search_h */