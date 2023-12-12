#ifndef SPHERE
#define SPHERE

#include "Vector3D.h"
#include <GL/glut.h>

#define pi (2*acos(0.0))

class Sphere
{
public:
    Vector3D *position = new Vector3D();
    Vector3D *up = new Vector3D();
    Vector3D *center = new Vector3D();
    Vector3D *look = new Vector3D();
    Vector3D *right = new Vector3D();  
    Vector3D *forward = new Vector3D();

    double angular_speed = .5, linear_speed = .5, radius = 1;
    int stacks = 40, slices = 40;
    double angular_rotation = 0;

    Sphere(Vector3D *position, Vector3D *up, Vector3D *center){
        this->position = position;
        position->z += radius;  // move the sphere center radius amount along z-axis

        this->center = center;  // where is our focus ,[center = position - look]
        this->up = up;

        *this->look = *this->position - *this->center; // front of the camera
        *this->right = ((*this->up) * (*this->look)).normalize_vector(); // up, look and right are perpendicular to each other
        *this->forward = (*this->right);        // current direction of the camera
    }

    ~Sphere(){
        delete position;
        delete up;
        delete center;
        delete look;
        delete right;
        delete forward;
    }

    void sphere_part(){
        Vector3D points[100][100];
        int i,j;
        double h,r;
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices/4;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            // glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
            for(j=0;j<slices/4;j++)
            {
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                }glEnd();
            }
        }
    }

    void drawSphereUsingParts(){
        glColor3f(0.5, 0, 0);
        glPushMatrix();
        {
            sphere_part ();
        }
        glPopMatrix();

        // Draw second sphere part
        glColor3f(0, 0.5, 0);
        glPushMatrix();
        {
            glRotatef(90,0,0,1);
            sphere_part ();
        }
        glPopMatrix();

        glColor3f(0, 0, 0.5);
        glPushMatrix();
        {
            glRotatef(180,0,0,1);
            sphere_part ( );
        }
        glPopMatrix();

        glColor3f(0, 0.5, 0.5);
        glPushMatrix();
        {
            glRotatef(270,0,0,1);
            sphere_part ( );
        }
        glPopMatrix();

        glColor3f(0.5, 0, 0.5);
        glPushMatrix();
        {
            glRotatef(90,0,1,0);
            sphere_part ( );
        }
        glPopMatrix();

        glColor3f(0.5, 0.2, 0.4);
        glPushMatrix();
        {
            glRotatef(90,0,0,1);
            glRotatef(90,0,1,0);
            sphere_part ( );
        }
        glPopMatrix();

        glColor3f(0.5, 0.5, 0);
        glPushMatrix();
        {
            glRotatef(180,0,0,1);
            glRotatef(90,0,1,0);
            sphere_part ( );
        }
        glPopMatrix();

        glColor3f(0, 1, 0.5);
        glPushMatrix();
        {
            glRotatef(270,0,0,1);
            glRotatef(90,0,1,0);
            sphere_part ( );
        }
        glPopMatrix();

        
    }

    void drawSphereUsingStacks(){

        // for rolling purpose
        glRotatef(angular_rotation, up->x, up->y, up->z);

        // stack: no of divisoin along y-axis(vertical axis)
        // slice: no of division along x-axis(horizontal axis)
        for (int i = 0; i <= stacks; ++i) {
            // draw a series of quadrilateral strips to form the sphere
            glBegin(GL_QUAD_STRIP);

            for (int j = 0; j <= slices; ++j) {
                float phi = i * pi / stacks;
                float theta = j * 2 * pi / slices;

                if(i<stacks/2){
                    (j%2) ? glColor3f(1,0,0): glColor3f(0,1,0);
                }
                else{
                    (j%2) ? glColor3f(0,1,0): glColor3f(1,0,0);
                }

                float x = radius * sin(phi) * cos(theta);
                float y = radius * sin(phi) * sin(theta);
                float z = radius * cos(phi);

                glVertex3f(x, y, z);

                phi = (i + 1) * M_PI / stacks;
                x = radius * sin(phi) * cos(theta);
                y = radius * sin(phi) * sin(theta);
                z = radius * cos(phi);

                glVertex3f(x, y, z);
            }
        }

        glEnd();

    }

    void draw_arrow(){
        glPushMatrix();
        // Draw a cylinder as the arrow shaft
        gluCylinder(gluNewQuadric(),.1, .1, 2.0, slices, slices);

        // Move the position to the top of the cylinder for the arrowhead
        glTranslatef(0.0, 0.0, 2.0);

        // Draw a cone as the arrowhead
        glutSolidCone(0.15, 0.25, 20, 20);

        glPopMatrix();
    }

    void drawAll(){ 
        glPushMatrix();{ 
            // move the origin sphere radius amount along z-axis 
            glTranslatef(position->x, position->y, position->z);
            drawSphereUsingStacks();
        }glPopMatrix();
        
        
        /////// Draw the upWard arrow ///////////////////////
        glPushMatrix();{
            glTranslatef(position->x, position->y, position->z+radius);
            glColor3f(0,1,0);
            draw_arrow();
        }glPopMatrix();

        /////// Draw the Direction arrow ///////////////////////
        glPushMatrix();{
            // Translate to the back of the sphere
            glTranslatef(position->x, position->y,position->z);
            // Calculate the rotation angles to align the arrow with the forward vector
            double angle = acos(forward->z);  // Angle between the forward vector and the z-axis
            double axisX = -forward->y;        // x-component of the cross product of forward and z-axis
            double axisY = forward->x;         // y-component of the cross product of forward and z-axis
            glRotatef(angle * 180.0 / M_PI, axisX, axisY, 0.0);
            glColor3f(0,0,1);
            draw_arrow();
        }glPopMatrix();
    }

    // press 'i' to move forward
    void move_forward(){
        cout<<"move_forward called\n";
        // ensure in the range [0, 360]
        while(angular_rotation > 360) angular_rotation -= 360;

        // convert to degree
        angular_rotation += ((linear_speed / radius) * (180.0/M_PI));


        *this->position = *this->position + ((*this->forward) * linear_speed); 
        

        // update look as the position is changed
        *this->look = (*this->right * (*this->up)).normalize_vector(); 


        
        check();        // check if the end of the block is reached
    }

    // press 'k' to move backward
    void move_backward(){
        cout<<"move backward called\n";
        // ensure in the range [0, 360]
        while(angular_rotation < -360) angular_rotation += 360;

        // convert to degree
        angular_rotation += ((linear_speed / radius) * (180.0/pi));

        // update position
        *this->position = *this->position - (*this->forward * linear_speed); 

        // update look as the position is changed
        *this->look = (*this->right * (*this->up)).normalize_vector();        // up, look and right are perpendicular to each other

        check();        // check if the end of the block is reached 
    }

    // key j pressed
    void rotate_arrow_anticlockwise(){
        printf("rotate_arrow_anticlockwise called\n");
        Vector3D temp = (*this->forward * *this->look).normalize_vector();  
        temp = temp * -tan(angular_speed);  // scale temp by tan(angular_speed)
        *this->forward = (*this->forward + &temp).normalize_vector();    // update forward
    }

    // key l pressed
    void rotate_arrow_clockwise(){
        printf("rotate_arrow_clockwise called\n");
        Vector3D temp = (*this->forward * *this->look).normalize_vector();   // perpendicular to forward and up
        temp = temp * tan(angular_speed);  // scale temp by tan(angular_speed)
        *this->forward = (*this->forward + &temp).normalize_vector();    // update forward

    }

    void check(){
        // since each red tile has size = 0.1 and there are 51 tiles
        double border = 5.0;

        if(position->x + radius >= (border-.2)) new_direction("right");
        else if(position->y + radius >= (border)) new_direction("top");
        else if(position->x - radius <= -(border-.2)) new_direction("left");
        else if(position->y - radius <= -(border)) new_direction("bottom");        
    }

    void new_direction(string s){
        if(s == "top" || s == "bottom"){
            forward->y = -forward->y;
        }
        if(s == "left" || s == "right"){
            forward->x = -forward->x;
        }
    }

    void print_all(){
        cout<<"position: "; position->print_vector();
        cout<<"up: "; up->print_vector();
        cout<<"look: "; look->print_vector();
        cout<<"right: "; right->print_vector();
        cout<<"forward: "; forward->print_vector();
        cout<<endl;
    }


};


#endif // SPHERE  