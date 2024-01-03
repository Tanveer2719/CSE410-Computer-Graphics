#ifndef HOMOGENOUSVECTOR
#define HOMOGENOUSVECTOR

#include<bits/stdc++.h>
using namespace std;

class HVector
{
public:
    
    double x, y, z, w;

    HVector(){}

    HVector(double a, double b, double c){
        x = a;
        y = b;
        z = c;
        w = 1;
    }

    HVector(double a, double b, double c, double d){
        x = a;
        y = b;
        z = c;
        w = d;
    }
    
    HVector(HVector *v){
        x = v->x, y = v->y, z = v->z, w = v->w;
    }

    // operator overloading
    HVector operator +(const HVector &b) {
        return HVector(x+b.x, y+b.y, z+b.z);
    }
    HVector operator -(const HVector &b) {
        return HVector((x-b.x) , (y-b.y) , (z-b.z) );
    }
    HVector operator *(double b) { return HVector(x*b, y*b, z*b);}
    HVector operator /(double b) {return HVector(x/b, y/b, z/b);}
    // cross multiplication
    HVector operator *(const HVector &b) {return HVector(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);}

    // for input
    friend istream& operator >>(istream& is, HVector& vec){
        is>> vec.x>> vec.y>> vec.z;
        vec.w = 1;
        return is;
    }

    // for console output
    friend ostream& operator << (ostream& os, HVector& vec){
        os << vec.x <<" "<< vec.y<<" "<< vec.z;
        return os;
    }

    // for file output
    friend ofstream& operator << (ofstream& os, HVector& vec){
        os << fixed << setprecision(7) << vec.x <<" "<< vec.y<<" "<< vec.z;
        return os;
    }
   
    HVector* normalize_vector(){
        // the normaliztion does not affect the w 
        // the normaliztion done on other points
		double magnitude = sqrt(x*x + y*y + z*z);
        if (! magnitude) return this;
		x /= magnitude;
		y /= magnitude;
		z /= magnitude;
        w = 1;
        return this;
	}

    HVector* scale_down(){
        x /= w;
        y /= w;
        z /= w;
        w /= w;
        return this;
    }

	void print_vector(){
		printf("x: %lf, y: %lf, z: %lf, w: %lf\n", x, y, z, w);
	}

    double dot(HVector *b){
        return x*b->x + y*b->y + z*b->z;
    }
};

#endif