//
//  Interpolate.h
//  pairproduce
//
//  Created by s610v on 7/29/16.
//  Copyright © 2016 s610v. All rights reserved.
//

#ifndef Interpolate_h
#define Interpolate_h

#include <stdio.h>
#include "IntFunction.h"
#include "search.h"

class Interpolate{
public:
    void interpolation(double*, double, Search, IntFunction);
    void polint(double*, double*, int, double, double*, double*);
};

#endif /* Interpolate_h */
