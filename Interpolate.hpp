//
//  Interpolate.hpp
//  pairproduce
//
//  Created by s610v on 7/19/16.
//  Copyright © 2016 s610v. All rights reserved.
//

#ifndef Interpolate_hpp
#define Interpolate_hpp

#include <stdio.h>
#include "IntFunction.h"
#include "search.h"

class Interpolate{
    float w2I;
public:
    Interpolate();
    Interpolate(float);
    void interpolation(float*, Search, IntFunction);
    void polint(float*, float*, int, float, float*, float*);
};

#endif /* Interpolate_hpp */
