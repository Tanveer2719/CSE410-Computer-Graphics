#ifndef MATRIX_H
#define MATRIX_H

#include "HomogenousVector.h"


class Matrix{
public:
    vector<vector<double>> mat;
    int dimension;

    Matrix(){
        mat.resize(4,vector<double>(4,0));
        dimension = 4;
        create_I();
    }
    
    Matrix(int d){
        mat.resize(d,vector<double>(d,0));
        dimension = d;
        create_I();
    }

    Matrix* get_matrix(){
        return this;
    }
    
    // Indentity Matrix
    void create_I(){
        for(int i=0;i<dimension;i++){
            for(int j=0;j<dimension;j++){
                if(i==j) mat[i][j] = 1;
                else mat[i][j] = 0;
            }
        }
    }

    // translation Matrix
    void create_T(HVector p){
        create_I();
        mat[0][dimension-1] = p.x;
        mat[1][dimension-1] = p.y;
        mat[2][dimension-1] = p.z;
    }

    // scaling matrix
    void create_S(HVector p){
        create_I();
        mat[0][0] = p.x;
        mat[1][1] = p.y;
        mat[2][2] = p.z;
    }

    // rodrigues
    HVector rodrigues(HVector &x, HVector &a, double angle){
        // angle in radian
        double ang = angle * M_PI / 180.0;
        return x*cos(ang) + a*(a.dot(&x))*(1-cos(ang)) + (a*x)*sin(ang);
    }

    // rotation matrix
    void create_R(HVector& a, double angle){
        create_I();

        auto i = HVector(1,0,0);
        auto j = HVector(0,1,0);
        auto k = HVector(0,0,1);
        
        a = a.normalize_vector();
        auto c1 = rodrigues(i, a, angle); // i
        auto c2 = rodrigues(j, a, angle); // j
        auto c3 = rodrigues(k, a, angle); // k
  
        mat[0][0] = c1.x, mat[0][1] = c2.x, mat[0][2] = c3.x,
        mat[1][0] = c1.y, mat[1][1] = c2.y, mat[1][2] = c3.y,
        mat[2][0] = c1.z, mat[2][1] = c2.z, mat[2][2] = c3.z;

    }

    // matrix vector multiplication
    // 1x4 * 4x4 = 1x4
    HVector operator *(HVector& b){
        double res[4];
        double b_vec[4] = {b.x, b.y, b.z, b.w};
        for(int i=0;i<dimension;i++){
            res[i] = 0;
            for(int j=0;j<dimension;j++){
                res[i] += mat[i][j]*b_vec[j];
            }
        }
        HVector ret(res[0],res[1],res[2],res[3]);
        ret.scale_down();
        return ret;
    }

    //matrix multiplication
    //4x4 * 4x4 = 4x4
    Matrix operator *(Matrix &b){
        Matrix ret (dimension);
        for(int i=0;i<dimension;i++){
            for(int j=0;j<dimension;j++){
                ret.mat[i][j] = 0;
                for(int k=0;k<dimension;k++){
                    ret.mat[i][j] += mat[i][k]*b.mat[k][j];
                }
            }
        }
        return ret;
    }


    void create_V(HVector &eye, HVector &look, HVector &up){
        create_I();

        auto l = (look - eye);
        l.normalize_vector();
        auto r = (l * up);
        r.normalize_vector();
        auto u = (r * l);
        u.normalize_vector();

        // translation T to move the eye/camera to origin.
        Matrix translate = Matrix(4);
        translate.mat[0][dimension-1] = -eye.x;
        translate.mat[1][dimension-1] = -eye.y;
        translate.mat[2][dimension-1] = -eye.z;

        // the l aligns with the -Z axis, r with X axis, and u with Y axis
        Matrix rotate = Matrix(4);
        rotate.mat[0][0] = r.x, rotate.mat[0][1] = r.y, rotate.mat[0][2] = r.z,
        rotate.mat[1][0] = u.x, rotate.mat[1][1] = u.y, rotate.mat[1][2] = u.z,
        rotate.mat[2][0] = -l.x, rotate.mat[2][1] = -l.y, rotate.mat[2][2] = -l.z;

        // translate.print_matrix();
        // cout<<endl; 
        // rotate.print_matrix();

        // view matrix
        mat = (rotate * translate).mat;
    }

    void create_P(double fovY, double aspectRatio, double near, double far){
        create_I();
        double fovX = fovY * aspectRatio;
        double t = near * tan((fovY/2) * M_PI / 180.0);
        double r = near * tan((fovX/2) * M_PI / 180.0);

        mat[0][0] = near/r,
        mat[1][1] = near/t,
        mat[2][2] = -(far+near)/(far-near),
        mat[2][3] = -(2*far*near)/(far-near),
        mat[3][2] = -1, mat[3][3] = 0;
    }

    void print_matrix(){
        for(int i=0;i<dimension;i++){
            for(int j=0;j<dimension;j++){
                cout<<mat[i][j]<<" ";
            }
            cout<<endl;
        }
    }
};

#endif