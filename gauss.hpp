//
//  gauss.hpp
//  pairproduce
//
//  Created by Vijay Shah on 7/25/16.
//  Copyright © 2016 Vijay Shah. All rights reserved.
//

#ifndef gauss_hpp
#define gauss_hpp

#include <stdio.h>
#define FUNC(x) ((*func)(x))

template <class T>
float qgaus(const T& obj, float (T::*func) (float) const, float a, float b){
    static const float x[] = {0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717};
    static const float w[] = {0.2955242247147529, 0.2692667193099963, 0.2190863625159821, 0.1494513491505806, 0.0666713443086881};
    float xm = 0.5*(b+a);
    float xr = 0.5*(b-a);
    float s = 0;
    for (int j=0; j<5; j++) {
        float dx = xr*x[j];
        s += w[j]*( (obj.*func)(xm+dx)+(obj.*func)(xm-dx));
    }
    return s *= xr;
}

#endif /* gauss_hpp */
