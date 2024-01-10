#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "HomogenousVector.h"

#define inf 1e9
#define epsilon 1e-9

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
    double a, b, c, d;  // for determining the plane equation

    Triangle(){
        setcolors();
    }
    Triangle(HVector a, HVector b, HVector c){
        vectors[0] = a;
        vectors[1] = b;
        vectors[2] = c;
        setcolors();
        determine_plane_equation();
    }

    void setcolors(){
        for(int i = 0; i<3; i++){
            colors[i] = ran_gen() % 255+1;
        }   
    }

    void determine_plane_equation(){
        HVector v1, v2;
        v1.x = vectors[1].x - vectors[0].x;
        v1.y = vectors[1].y - vectors[0].y;
        v1.z = vectors[1].z - vectors[0].z;

        v2.x = vectors[2].x - vectors[0].x;
        v2.y = vectors[2].y - vectors[0].y;
        v2.z = vectors[2].z - vectors[0].z;

        HVector normal = v1*v2;
        // cout<<"normal: "<<normal<<endl;
        normal.normalize_vector();
        
        // cout<<"v1: "<<v1<<endl;
        // cout<<"v2: "<<v2<<endl;
        // cout<<"normal :"<<normal<<endl;
    

        this->a = normal.x;
        this->b = normal.y;
        this->c = normal.z;
        this->d = -(a*vectors[0].x + b*vectors[0].y + c*vectors[0].z);

    }

    void sort(){
        if(vectors[0].y > vectors[1].y) swap(vectors[0], vectors[1]);
        if(vectors[0].y > vectors[2].y) swap(vectors[0], vectors[2]);
        if(vectors[1].y > vectors[2].y) swap(vectors[2], vectors[1]);

    }

    pair<double, double> find_left_right_intersection(double y, double left_x, double right_x, double dx, double dy){
        pair<double, double> intersections(99999, -99999); // minimum intersection, maximum intersection

        for(int i = 0; i<3; i++){
            int j = (i+1)%3;

            double x1 = vectors[i].x, x2 = vectors[j].x;
            double y1 = vectors[i].y, y2 = vectors[j].y;

            if(y1 > y2){
                swap(x1, x2); 
                swap(y1, y2);
            }

            if(y >= y1 && y<= y2 && abs(y2-y1) > epsilon){
                double x = x1 + (x2-x1)*(y-y1)/(y2-y1);
                intersections.first = min(intersections.first, x);
                intersections.second = max(intersections.second, x);

            }
        
        }
        

        intersections.first = round(max(intersections.first/dx, left_x/dx));
        intersections.second = round(min(intersections.second/dx, right_x/dx));
        return intersections;
    }

};

#endif