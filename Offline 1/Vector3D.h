#ifndef VECTOR3D
#define VECTOR3D

#include<bits/stdc++.h>
using namespace std;

class Vector3D
{
public:
    
    double x, y, z;

    Vector3D(){}

    Vector3D(double a, double b, double c){
        x = a;
        y = b;
        z = c;
    }
    
    Vector3D(Vector3D *v){
        x = v->x, y = v->y, z = v->z;
    }

    // operator overloading
    Vector3D operator +(const Vector3D &b) {return Vector3D(x+b.x, y+b.y, z+b.z);}
    Vector3D operator -(const Vector3D &b) {return Vector3D(x-b.x, y-b.y, z-b.z);}
    Vector3D operator *(double b) {return Vector3D(x*b, y*b, z*b);}
    Vector3D operator /(double b) {return Vector3D(x/b, y/b, z/b);}

    // cross multiplication
    Vector3D operator *(const Vector3D &b) {return Vector3D(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);}

	Vector3D* normalize_vector(){
		double magnitude = sqrt(x*x + y*y + z*z);
        if (! magnitude) return this;
		x /= magnitude;
		y /= magnitude;
		z /= magnitude;
        return this;
	}

	void print_vector(){
		printf("x: %lf, y: %lf, z: %lf\n", x, y, z);
	}

    Vector3D* dot(Vector3D *b){
        return new Vector3D(x*b->x, y*b->y, z*b->z);
    }
};

#endif