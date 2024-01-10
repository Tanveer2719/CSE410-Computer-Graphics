#include "HomogenousVector.h"
#include "Matrix.h"
#include "Triangle.h"
#include "bitmap_image.hpp"

using namespace std;

int get_left_intersection(const Triangle& t, double y, double left_x, double dx, int triangle_no){
    double min_x = 99999;

    for(int i = 0; i<3; i++){
        int j = (i+1)%3;
        double x1 = t.vectors[i].x;
        double y1 = t.vectors[i].y;
        double x2 = t.vectors[j].x;
        double y2 = t.vectors[j].y;

        if(y1 > y2){
            swap(x1, x2);
            swap(y1, y2);  
        }

        if(y >= y1 and y<= y2 and (y2-y1)>epsilon){
            double x = x1 + (y-y1)*(x2-x1)/(y2-y1);
            min_x = min(min_x, x);
        }
    }
    
    int xa = round(max(min_x/dx, left_x/dx));
    return xa;
}

int get_right_intersection(const Triangle& t, double y, double right_x, double dx, int triangle_no){
    double max_x = -99999;

    for(int i = 0; i<3; i++){
        int j = (i+1)%3;
        double x1 = t.vectors[i].x;
        double y1 = t.vectors[i].y;
        double x2 = t.vectors[j].x;
        double y2 = t.vectors[j].y;

        if(y1 > y2){
            swap(x1, x2);
            swap(y1, y2);  
        }

        if(y >= y1 and y<= y2 and (y2-y1)>epsilon){
            double x = x1 + (y-y1)*(x2-x1)/(y2-y1);
            max_x = max(max_x, x);
        }
    }
 
    int xb = round(min(max_x/dx, right_x/dx));
    return xb;
}

int main(){

    HVector eye, look, up;
    double fov, asp_r, near, far;
    int n_triangles = 0;
    


    ifstream input_file("5/scene.txt");
    ofstream out_file("stage1.txt");

    if(input_file && out_file){
        cout<<"Both files opened successfully"<<endl;   
    }else{
        cout<<"error in file opening"<<endl;
        return 0;
    }   

    input_file >> eye >> look >> up;
    input_file >> fov >> asp_r >> near >> far;

    
    /*
        @brief : stage 1
    */

    //stack
    stack<Matrix> stack;
    Matrix c_matrix (4);
    stack.push(c_matrix);

    while (true){
        string command;
        input_file >> command;

        if(command == "triangle"){
            HVector v1, v2, v3;
            input_file >> v1 >> v2 >> v3;


            v1 = (c_matrix) * (v1);
            v2 = (c_matrix) * (v2);
            v3 = (c_matrix) * (v3);


            Triangle t = Triangle(v1, v2, v3);
            n_triangles++;

            out_file << v1 << endl;
            out_file << v2 << endl;
            out_file << v3 << endl; 
            out_file << endl;
        }
        
        else if(command == "end"){
            break;
        }
        
        else if(command == "translate"){
            double tx, ty, tz;
            input_file >> tx >> ty >> tz;

            Matrix t_matrix = Matrix(4);
            t_matrix.create_T(HVector(tx, ty, tz));
            c_matrix = (c_matrix) * (t_matrix);
        }
        
        else if(command == "scale"){
            double sx, sy, sz;
            input_file >> sx >> sy >> sz;

            Matrix s_matrix = Matrix(4);
            s_matrix.create_S(HVector(sx, sy, sz));
            c_matrix = (c_matrix) * (s_matrix);
        }

        else if(command == "rotate"){
            double angle, ax, ay, az;
            input_file >> angle >> ax >> ay >> az;  
            HVector vec  = HVector(ax, ay,az);

            Matrix r_matrix = Matrix(4);
            r_matrix.create_R(vec, angle);
            c_matrix = (c_matrix) * (r_matrix);
        }
    
        else if(command == "push"){
            stack.push(c_matrix);
        }
        
        else if(command == "pop"){
            if(stack.empty()){
                cout << "Stack is empty" << endl;
                c_matrix.create_I(); 
                stack.push(c_matrix);   
            }
            c_matrix = stack.top();
            stack.pop(); 
        }
    }

    input_file.close();
    out_file.close();

    /*
        @brief : stage 2
    */

    input_file.open("stage1.txt");
    out_file.open("stage2.txt");
    if(out_file){
        cout<<"stage2.txt opened successfully"<<endl;
    }
    else{
        cout<<"error in stage2.txt opening"<<endl;
        return 0;
    }

    Matrix view_T (4);
    view_T.create_V(eye, look, up);

    for(int i = 0; i<n_triangles; i++){
        HVector v1, v2, v3;
        input_file >> v1 >> v2 >> v3;
        Triangle t = Triangle(v1, v2, v3);
        t.vectors[0] = (view_T) * (t.vectors[0]);
        t.vectors[1] = (view_T) * (t.vectors[1]);
        t.vectors[2] = (view_T) * (t.vectors[2]);

        out_file << t.vectors[0] << endl;
        out_file << t.vectors[1] << endl;
        out_file << t.vectors[2] << endl; 
        out_file << endl;
    }
    input_file.close();
    out_file.close();

    /*
        @brief : stage 3
    */  
    out_file.open("stage3.txt");
    if(out_file){
        cout<<"stage3.txt opened successfully"<<endl;
    }
    else{
        cout<<"error in stage3.txt opening"<<endl;
        return 0;
    }

    Matrix view_P (4);
    view_P.create_P(fov, asp_r, near, far);

    input_file.open("stage2.txt");
    for(int i = 0; i<n_triangles; i++){
        HVector v1, v2, v3;
        input_file >> v1 >> v2 >> v3;
        Triangle t = Triangle(v1, v2, v3);
        t.vectors[0] = (view_P) * (t.vectors[0]);
        t.vectors[1] = (view_P) * (t.vectors[1]);
        t.vectors[2] = (view_P) * (t.vectors[2]);

        out_file << t.vectors[0] << endl;
        out_file << t.vectors[1] << endl;
        out_file << t.vectors[2] << endl; 
        out_file << endl;
    }

    input_file.close();
    out_file.close();

    /*
        @brief : stage 4
    */ 

    double s_width, s_height;
    input_file.open("config.txt");
    out_file.open("z_buffer.txt");

    // read
    input_file >> s_width >> s_height;
    input_file.close();
    
    input_file.open("stage3.txt");

    // right-left-top-bottom limits
    double rl = 1, ll = -1, tl = 1, bl = -1;
    double z_max = 1.0, z_min = -1.0;
    

    // initialize z-buffer and frame buffer    
    // dx*dy =  area of a cell of the screen(pixel)
    double dx = (rl - ll) / s_width;
    double dy = (tl - bl) / s_height;
    
    double top_y = tl - dy/2;
    double bottom_y = bl + dy/2;
    double left_x = ll + dx/2;
    double right_x = rl - dx/2; 

    // cout<<"top_y: "<<top_y<<" bottom_y: "<< bottom_y<<" left_x: "<<left_x<<" right_x: "<<right_x<<endl;


    // initialize z-buffer
    vector<vector<double> > z_buffer(s_height, vector<double>(s_width, z_max));

    for(int i = 0; i<s_height; i++){
        for(int j = 0; j<s_width; j++){
            z_buffer[i][j] = z_max;
        }
    }

    // initialize image
    bitmap_image image(s_width, s_height);
    
    // set color to black
    image.set_all_channels(0,0,0);

    // read triangles
    for(int k = 0; k<n_triangles; k++){
        HVector v1, v2, v3;
        input_file >> v1 >> v2 >> v3;

        Triangle t = Triangle(v1, v2, v3);
        t.determine_plane_equation();
        t.sort();
    
        // cout<<"*************** triangle "<<k+1 <<"***************"<<endl;
        // cout<<v1<<endl;
        // cout<<v2<<endl;
        // cout<<v3<<endl;
        // cout<<"*************** triangle ***************"<<endl;


        // clip y_axis
        int  min_y = floor(max(t.vectors[0].y, bottom_y)/dy);
        int  max_y = ceil(min(t.vectors[2].y, top_y)/dy);

        
        for(int y_val = max_y; y_val>= min_y; y_val--){
 
            double y = y_val*dy;
            
            int xa = get_left_intersection(t, y, left_x, dx, k); 
            int xb = get_right_intersection(t, y, right_x, dx, k);
            
            for(int x_val = xa; x_val<= xb; x_val++){
                double x = x_val*dx;

                int i = top_y/dy - y_val;
                int j = x_val - left_x/dx;

                //calculate z-values
                double zp = - (t.a*x + t.b*y + t.d) / t.c;


                if(zp < z_buffer[i][j] && zp > z_min){
                    z_buffer[i][j] = zp;
                    image.set_pixel(j, i, t.colors[0], t.colors[1], t.colors[2]);
                }
            }
        }

    } 
    
    for(int i = 0; i<s_height; i++){
        for(int j = 0; j<s_width; j++){
            if(z_buffer[i][j] < z_max){
                out_file<<fixed<<setprecision(6)<<z_buffer[i][j]<<"\t";
            }
        }
        out_file<<endl;
    }

    out_file.close();
    input_file.close();

    image.save_image("out.bmp");
    z_buffer.clear();
    z_buffer.shrink_to_fit();

    return 0;
}