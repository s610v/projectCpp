//
//  gauss.hpp
//  pairproduce
//
//  Created by s610v on 7/25/16.
//  Copyright © 2016 s610v. All rights reserved.
//

#ifndef gauss_hpp
#define gauss_hpp

#include <stdio.h>
#define FUNC(x) ((*func)(x))

template <class T>
double qgaus(const T& obj, double (T::*func) (double) const, double a, double b){
    static const double x[] = {0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717};
    static const double w[] = {0.2955242247147529, 0.2692667193099963, 0.2190863625159821, 0.1494513491505806, 0.0666713443086881};
    double xm = 0.5*(b+a);
    double xr = 0.5*(b-a);
    double s = 0;
    for (int j=0; j<5; j++) {
        double dx = xr*x[j];
        s += w[j]*( (obj.*func)(xm+dx)+(obj.*func)(xm-dx));
    }
    return s *= xr;
}

#endif /* gauss_hpp */
