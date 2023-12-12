#ifndef CAMERA
#define CAMERA

#include "Vector3D.h"

class Camera
{
public:
    Vector3D *position = new Vector3D();    // where the camera is
    Vector3D *look = new Vector3D();        // where the camera is facing
    Vector3D *up = new Vector3D();          // the up vector
    Vector3D *right = new Vector3D();       // the right vector
    Vector3D *center  = new Vector3D();     // currently where the camera is looking at

    double linear_speed = .3, angular_speed = .05;

    Camera(Vector3D *position, Vector3D *center, Vector3D *up){
        this->position = position;
        this->center = center;
        this->up = up;

        *this->look = *this->position - *this->center;
        *this->right = (*this->up) * (*this->look);         // up, look and right are perpendicular to each other
        *this->right->normalize_vector();
    }

    // assume that the direction is normalized
    // need to update the look and the position
    void move_camera(Vector3D * direction) { 
        *this->position = (*this->position + *direction);
        *this->center = (*this->center + *direction);
    }

    //The logic involves \
    subtracting the current position (pos) \
    from the target, normalizing the result, \
    multiplying it by the magnitude of the \
    original direction, and then adding the \
    position back

    void rotate_camera(Vector3D* direction){    
        Vector3D temp = *this->center - *this->position;
        double magnitude = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
        
        temp = temp.normalize_vector();

        temp = temp * magnitude;
        *this->center = temp + *direction + *this->position;
    }

    // up arrow key
    void move_forward(){
        Vector3D direction = *this->look * -linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }
    
    // down arrow key
    void move_backward(){
        Vector3D direction = *this->look * linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }

    // right arrow key
    void move_right(){
        Vector3D direction = *this->right * (-linear_speed);
        direction.normalize_vector();
        move_camera(&direction);
    }

    // left arrow key
    void move_left(){
        Vector3D direction = *this->right * (linear_speed);
        direction.normalize_vector();
        move_camera(&direction);
    }

    // page up key
    void move_up(){
        Vector3D direction = *this->up * linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }

    // page down key
    void move_down(){
        Vector3D direction = *this->up * -linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }

    // 3 key
    void look_up(){
        Vector3D direction = *this->up * (tan(angular_speed)); 
        rotate_camera(&direction);
    }

    // 4 key
    void look_down(){
        Vector3D direction = *this->up * (-tan(angular_speed)); 
        rotate_camera(&direction);
    }


    // 1 key
    void look_left(){
        Vector3D direction = *this->right * (-tan(angular_speed)); 
        rotate_camera(&direction);
    }

    // 2 key
    void look_right(){
        Vector3D direction = *this->right * (tan(angular_speed)); 
        rotate_camera(&direction);  
    }

    // 6 key
    void tilt_clockwise(){
        *this->up = ((*this->right) 
                    * (-tan(angular_speed)) 
                    + (*this->up))
                    .normalize_vector();

        *this->right = ( (*this->up) 
                        * (*this->look))
                        .normalize_vector();

    }
    // 5 key
    void tilt_anticlokwise(){

        *this->up = ((*this->right) 
                    * tan(angular_speed) 
                    + (*this->up))
                    .normalize_vector();

        *this->right = ((*this->up) 
                        * (*this->look))
                        .normalize_vector();
         
    }

    // w key
    void up_reference(){
        *this->position = ((*this->up) 
                    * (linear_speed))
                    + (*this->position)
                    ;
    }

    // s key
    void down_reference(){
        *this->position = ((*this->up) 
                    * (-linear_speed))
                    + (*this->position)
                   ;
    }





    ~Camera(){
        delete position;
        delete look;
        delete up;
        delete right;
        delete center;
    }
    

};




#endif
