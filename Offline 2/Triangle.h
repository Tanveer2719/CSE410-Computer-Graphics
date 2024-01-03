#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "HomogenousVector.h"

static unsigned long int g_seed = 1;
inline int ran_gen()
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

class Triangle{
public:
    HVector vectors[3];
    double colors[3];

    Triangle(){
        setcolors();
    }
    Triangle(HVector a, HVector b, HVector c){
        vectors[0] = a;
        vectors[1] = b;
        vectors[2] = c;
        setcolors();
    }
    void setcolors(){
        for(int i = 0; i<3; i++){
            colors[i] = ran_gen() % 255+1;
        }   
    }

};

#endif