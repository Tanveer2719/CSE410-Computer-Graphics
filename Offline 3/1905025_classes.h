#ifndef CLASSES
#define CLASSES

#include<bits/stdc++.h>
#include "bitmap_image.hpp"
#include<GL/glut.h>
#define pi (2*acos(0.0))
#define epsilon 1e-4
using namespace std;

class Vector3D 
{
public:
    
    double x, y, z, w;

    Vector3D (){}

    Vector3D (double a, double b, double c){
        x = a;
        y = b;
        z = c;
        w = 1;
    }

    Vector3D (double a, double b, double c, double d){
        x = a;
        y = b;
        z = c;
        w = d;
    }
    
    Vector3D (Vector3D  *v){
        x = v->x, y = v->y, z = v->z, w = v->w;
    }

    // operator overloading
    Vector3D  operator +(const Vector3D  &b) {
        return Vector3D (x+b.x, y+b.y, z+b.z);
    }
    Vector3D  operator -(const Vector3D  &b) {
        return Vector3D ((x-b.x) , (y-b.y) , (z-b.z) );
    }
    Vector3D  operator *(double b) {return Vector3D (x*b, y*b, z*b);}
    Vector3D  operator /(double b) {return Vector3D (x/b, y/b, z/b);}
    
    // cross multiplication
    Vector3D  operator *(const Vector3D  &b) {return Vector3D (y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x);}

    // negation
    Vector3D  operator -() {return Vector3D (-x, -y, -z);}
    
    // for input
    friend istream& operator >>(istream& is, Vector3D & vec){
        is>> vec.x>> vec.y>> vec.z;
        vec.w = 1;
        return is;
    }

    // for console output
    friend ostream& operator << (ostream& os, Vector3D & vec){
        os << vec.x <<" "<< vec.y<<" "<< vec.z;
        return os;
    }

    // for file output
    friend ofstream& operator << (ofstream& os, Vector3D & vec){
        os << fixed << setprecision(7) << vec.x <<" "<< vec.y<<" "<< vec.z;
        return os;
    }
   
    double magnitude(){
        return sqrt(x*x + y*y + z*z);
    }
   
    void normalize_vector(){
        // the normaliztion does not affect the w 
        // the normaliztion done on other points
		double m = magnitude();
        if (! m) return ;
		x /= m;
		y /= m;
		z /= m;
        w = 1;
	}

    void scale_down(){
        x /= w;
        y /= w;
        z /= w;
        w /= w;
    }

	void print_vector(){
		printf("x: %lf, y: %lf, z: %lf, w: %lf\n", x, y, z, w);
	}

    double dot(Vector3D  *b){
        return x*b->x + y*b->y + z*b->z;
    }

};

class Camera
{
public:
    Vector3D  *position = new Vector3D ();    // where the camera is
    Vector3D  *look = new Vector3D ();        // where the camera is facing
    Vector3D  *up = new Vector3D ();          // the up vector
    Vector3D  *right = new Vector3D ();       // the right vector
    Vector3D  *center  = new Vector3D ();     // currently where the camera is looking at

    double linear_speed = .3, angular_speed = .05;

    Camera(Vector3D * position, Vector3D * look, Vector3D * right){
        this->position = position;
        this->look = look;
        this->right = right;

        *this->center = *this->position + *this->look;
        *this->up = (*this->right) * (*this->look); 
    }

    // assume that the direction is normalized
    // need to update the look and the position
    void move_camera(Vector3D  * direction) { 
        *this->position = (*this->position + *direction);
        *this->center = (*this->center + *direction);
    }

    //The logic involves \
    subtracting the current position (pos) \
    from the target, normalizing the result, \
    multiplying it by the magnitude of the \
    original direction, and then adding the \
    position back

    void rotate_camera(Vector3D * direction){    
        Vector3D  temp = *this->center - *this->position;
        double magnitude = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
        
        temp.normalize_vector();

        temp = temp * magnitude;
        *this->center = temp + *direction + *this->position;
    }

    // up arrow key
    void move_forward(){
        Vector3D  direction = *this->look * linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }
    
    // down arrow key
    void move_backward(){
        Vector3D  direction = *this->look * -linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }

    // right arrow key
    void move_right(){
        Vector3D  direction = *this->right * (linear_speed);
        direction.normalize_vector();
        move_camera(&direction);
    }

    // left arrow key
    void move_left(){
        Vector3D  direction = *this->right * (-linear_speed);
        direction.normalize_vector();
        move_camera(&direction);
    }

    // page up key
    void move_up(){
        Vector3D  direction = *this->up * -linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }

    // page down key
    void move_down(){
        Vector3D  direction = *this->up * linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }

    // 3 key
    void look_up(){
        Vector3D  direction = *this->up * (-tan(angular_speed)); 
        rotate_camera(&direction);
    }

    // 4 key
    void look_down(){
        Vector3D  direction = *this->up * (+tan(angular_speed)); 
        rotate_camera(&direction);
    }


    // 1 key
    void look_left(){
        Vector3D  direction = *this->right * (tan(angular_speed)); 
        rotate_camera(&direction);
    }

    // 2 key
    void look_right(){
        Vector3D  direction = *this->right * (-tan(angular_speed)); 
        rotate_camera(&direction);  
    }

    // 6 key
    void tilt_clockwise(){
        *this->up = ((*this->right) 
                    * (tan(angular_speed)) 
                    + (*this->up));

        up->normalize_vector();

        *this->right = ( (*this->up) 
                        * (*this->look));
        right->normalize_vector(); 

    }
    // 5 key
    void tilt_anticlokwise(){

        *this->up = ((*this->right) 
                    * (-tan(angular_speed)) 
                    + (*this->up));
        up->normalize_vector();

        *this->right = ((*this->up) 
                        * (*this->look));
        right->normalize_vector();
         
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

// new classes


class Color {
public:
    double r, g, b;
    Color(){}
    Color(double r, double g, double b){
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color(const Color &c){
        r = c.r, g = c.g, b = c.b;
    }
    Color(const vector<double>&c){
        r=c[0], g=c[1], b=c[2];
    }

    friend ostream &operator <<(ostream& out, const Color &c){
        out<<c.r <<" "<<c.g<<" "<<c.b;
        return out;
    }

    // input from file
    friend istream &operator >>(istream &in, Color &c){
        in>>c.r>>c.g>>c.b;
        return in;
    }

    Color operator * (double d){
        return  Color(r*d, g*d, b*d);
    }
    Color operator + (const Color &c){
        return Color(r+c.r, g+c.g, b+c.b);
    }

};

class Ray{
public:
    Vector3D  origin, direction;
    Ray(){}
    Ray(Vector3D  origin, Vector3D  direction){
        this->origin = origin;
        direction.normalize_vector();
        this->direction = direction;
    }

    // // overload << for printing
    // friend ostream &operator << (ostream &out, const Ray &r){
    //     out<<"origin: "<<r.origin<<" direction: "<<r.direction<<"\n";
    //     return out;
    // }

};

class PointLight{
public:
    Vector3D  position;
    Color color;
    double power;       // intensity of the light source

    PointLight(){}
    void setColor(double c1, double c2, double c3){
        color = Color(c1, c2, c3);
    }
    void setColor(const vector<double> &color){
        this->color = Color(color);
    }
    void setColor(const Color &c){
        color = Color(c);
    }
    
    void setPower(double power){
        this->power = power;
    }
    void draw(){
        glPointSize(5);

        glBegin(GL_POINTS);
        glColor3f(color.r, color.g, color.b);
        glVertex3f(position.x, position.y, position.z);
        glEnd();
    }

    void printPointLight(){
        cout<<"PointLight\n";
        cout<<"position: "<<position<<endl;
        cout<<"color: "<<color<<endl;
    }

    // take input from file
    friend istream &operator >>(istream &in, PointLight &p){
        in>>p.position.x>>p.position.y>>p.position.z;
        in>>p.color.r>>p.color.g>>p.color.b;
        return in;
    }

    ~PointLight(){
        delete &position;
    }
};

class Spotlight{
public:
    PointLight pointLight;
    Vector3D  direction;
    double cutoff_angle;

    Spotlight(){}
    void setAngle(double angle){
        this->cutoff_angle = angle;
    }
    
    void setColor(double c1, double c2, double c3){
        pointLight.color = Color(c1, c2, c3);
    }
    
    void setColor(const vector<double> &color){
        pointLight.color = Color(color);
    }
    
    void draw(){
        glPointSize(15);
        glBegin(GL_POINTS);
        glColor3f(pointLight.color.r, pointLight.color.g, pointLight.color.b);
        glVertex3f(pointLight.position.x, pointLight.position.y, pointLight.position.z);
        glEnd();
    }

    // take input from file
    friend istream &operator >>(istream &in, Spotlight &s){
        in>>s.pointLight.position.x>>s.pointLight.position.y>>s.pointLight.position.z;
        in>>s.pointLight.color.r>>s.pointLight.color.g>>s.pointLight.color.b;
        in>>s.direction.x>>s.direction.y>>s.direction.z;
        in>>s.cutoff_angle;
        s.direction.normalize_vector();
        return in;
    }

    // print spotlight
    void printSpotlight(){
        cout<<"Spotlight\n";
        cout<<"position: "<<pointLight.position<<endl;
        cout<<"color: "<<pointLight.color<<endl;
        cout<<"direction: "<<direction<<endl;
        cout<<"cutoff_angle: "<<cutoff_angle<<endl;
    }

};

class Object;
extern vector<Object*>objects;
extern vector<PointLight*>  pointLights;
extern vector<Spotlight*> spotlights;

extern int recursionLevel;

//virtual Color getColor(Vector3D  p)
//virtual double objectRayIntersection(Ray &ray)=0;
//virtual void draw()=0;
//virtual Ray getNormal(Vector3D  point, Ray ray)=0;
//virtual double intersect(Ray ray, Color &color, int level)=0;
//virtual void print()=0;
//virtual ~Object(){}


class Object{

    bool is_obstructed(Vector3D  lightPosition, Ray incidentRay){
        double dist = (lightPosition - incidentRay.origin).magnitude();
        for(auto o : objects){
            double temp = o->objectRayIntersection(incidentRay);
            if(temp > 0 && temp + epsilon  < dist){
                return true;
            }
        }
        return false;
    }

    void update_color_components(Color &color, Color &source, Color &colorAtPoint, double lambert, double phong){
        color.r += source.r * lambert * colorAtPoint.r * coEfficients[1];
        color.r += source.r * pow(phong, shine) * colorAtPoint.r * coEfficients[2];

        color.g += source.g * lambert * colorAtPoint.g * coEfficients[1];
        color.g += source.g * pow(phong, shine) * colorAtPoint.g * coEfficients[2];

        color.b += source.b * lambert * colorAtPoint.b * coEfficients[1];
        color.b += source.b * pow(phong, shine) * colorAtPoint.b * coEfficients[2];
    }

    int find_nearest_object(Ray &ray, Color &color){
        int nearest = -1;
        double tMin = 1e9, t;
        for(int k=0;k<(int)objects.size();k++){
            t = objects[k]->intersect(ray, color, 0);
            if(t> 0 && t<tMin)
                tMin = t , nearest = k;
        }
        return nearest;
    }

public:
    Vector3D  reference_point; // should have x, y, z
    double height, width, length;
    Color color;
    vector<double>coEfficients; // ambient, diffuse, specular, reflection coEfficients
    int shine; // exponent term of specular component

    Object(){
        coEfficients.resize(4, 0.0);
    }
    
    void setColor(double c1, double c2, double c3){
        color = Color(c1, c2, c3);
    }
    
    void setColor(const vector<double> &colors){
        color = Color(colors);
        
    }
    
    void setColor(const Color &c){
        color = Color(c);
    }
    
    virtual Color getColor(Vector3D  p){
        return Color(this->color.r, 
                     this->color.g,
                     this->color.b);
    }

    void setShine(int exp){
        shine = exp;
    }
    
    void setcoEfficients(double ambient, double diffuse, double specular, double reflection){
        coEfficients[0] = ambient;
        coEfficients[1] = diffuse;
        coEfficients[2] = specular;
        coEfficients[3] = reflection;
    }
    
    void setcoEfficients(double coEfficients[]){
        for(int i = 0; i < 4; i++){
            this->coEfficients[i] = coEfficients[i];
        }
    }

    virtual double objectRayIntersection(Ray &ray)=0;
    
    virtual void draw()=0;
    
    virtual Ray getNormal(Vector3D  point, Ray ray)=0;

    virtual double intersect(Ray ray, Color &color, int level)
    {
        double t = objectRayIntersection(ray);
        if(t < 0) return -1;
        if(level == 0) return t;

        
        Vector3D  incidentPoint = ray.origin + ray.direction*t;
        Color colorAtIncidentPoint = getColor(incidentPoint);

        // ambient lighting
        color = colorAtIncidentPoint * coEfficients[0];

        for(int i = 0; i < pointLights.size(); i++){

            Vector3D  lightPosition = pointLights[i]->position;
            Vector3D  lightDirection = incidentPoint - lightPosition;
            lightDirection.normalize_vector();
            
            // incident ray
            Ray incidentRay = Ray(lightPosition, lightDirection);

            // calculate normal
            Ray normal = getNormal(incidentPoint,incidentRay);

            bool obscured = is_obstructed(lightPosition, incidentRay);
            
            if(!obscured){

                double inc_normal = incidentRay.direction.dot(&normal.direction);
                
                // lambert value
                double lambert = max(0.0, -inc_normal);
                
                // find reflected ray
                Ray reflectedRay = Ray(incidentPoint, incidentRay.direction - normal.direction*2*(inc_normal));
                double phong = max(0.0,-ray.direction.dot(&reflectedRay.direction));
                
                // update diffuse and specular components
                update_color_components(color, pointLights[i]->color, colorAtIncidentPoint, 
                                        lambert, phong);
            }
        }

        for(int i = 0; i < spotlights.size(); i++){
            Vector3D  lightPosition = spotlights[i]->pointLight.position;
            Vector3D  lightDirection = incidentPoint - lightPosition;
            lightDirection.normalize_vector();

            double temp = lightDirection.dot(&spotlights[i]->direction);
            double angle = acos(temp/(lightDirection.magnitude()*spotlights[i]->direction.magnitude())) * (180.0/pi);

            if(fabs(angle)<spotlights[i]->cutoff_angle){

                Ray incidentRay = Ray(lightPosition, lightDirection);
                Ray normal = getNormal(incidentPoint,incidentRay);
                
                bool obscured = is_obstructed(lightPosition, incidentRay);
                
                
                if(!obscured){
                    double inc_normal = incidentRay.direction.dot(&normal.direction);
                    Ray reflected_ray = Ray(incidentPoint, incidentRay.direction - normal.direction*2*(inc_normal));
                    
                    double phong = max(0.0,-ray.direction.dot(&reflected_ray.direction));
                    double lambert = max(0.0, -inc_normal);

                    update_color_components(color, spotlights[i]->pointLight.color, colorAtIncidentPoint, 
                                            lambert, phong);
                    
                }
            }
        }
 
        if(level < recursionLevel){

            // find normal at incidentPoint
            Ray normal = getNormal(incidentPoint,ray);

            // find reflected ray
            Ray reflectedRay = Ray(incidentPoint, ray.direction - normal.direction*2*(ray.direction.dot(&normal.direction)));

            // slightly shift the origin to avoid self intersection
            reflectedRay.origin = reflectedRay.origin + reflectedRay.direction*epsilon ;
            

            int nearest = find_nearest_object(reflectedRay, color);
            if(nearest != -1)
            {
                Color reflectionColor(0,0,0); // refelction color

                double t = objects[nearest]->intersect(reflectedRay,reflectionColor, level+1);

                color = color + reflectionColor*coEfficients[3];
            }
        }

        return t;
    }

/* 
    virtual double intersect(Ray ray, Color &color, int level)
    {
        double t = objectRayIntersection(ray);

        // cout<<"T "<<t<<endl;
        // return 0;

        if(t < 0) return -1;
        if(level == 0) return t;

        // find intersection point and it's color
        Vector3D intersectionPoint = ray.origin + ray.direction*t;
        Color colorAtIntersection = getColor(intersectionPoint);

        // update color with ambience (thing will become dimmer)
        color.r = colorAtIntersection.r * coEfficients[0];
        color.g = colorAtIntersection.g * coEfficients[0];
        color.b = colorAtIntersection.b * coEfficients[0];

        // cout<< " pointLights size " << pointLights.size() << endl;

        for(int i = 0; i < pointLights.size(); i++){

            Vector3D lightPosition = pointLights[i]->position;
            Vector3D lightDirection = intersectionPoint - lightPosition;
            lightDirection.normalize_vector();
            
            // cast incident ray, from light position to intersection point
            Ray lightRay = Ray(lightPosition, lightDirection);

            // calculate normal at intersectionPoint
            Ray normal = getNormal(intersectionPoint,lightRay);

    
            
             double t2 = (intersectionPoint - lightPosition).magnitude();
            if(t2 < 1e-5) continue;

            bool obscured = false;

            for(Object *obj : objects){
                double t3 = obj->objectRayIntersection(lightRay );
                if(t3 > 0 && t3 + 1e-5 < t2){
                    obscured = true;
                    break;
                }
            }

            if(!obscured){
                
                // lambert value
                double val = max(0.0, -lightRay.direction.dot(&normal.direction));
                
                // find reflected ray
                Ray reflection = Ray(intersectionPoint, lightRay.direction - normal.direction*2*(lightRay.direction.dot(&normal.direction)));
                double phong = max(0.0,-ray.direction.dot(&reflection.direction));
                
                // update diffuse and specular components
                // pointLights[i]->color works as the source intensity, Is here

                color.r += pointLights[i]->color.r * coEfficients[1] * val * colorAtIntersection.r;
                color.r += pointLights[i]->color.r * coEfficients[2] * pow(phong,shine) * colorAtIntersection.r;

                color.g += pointLights[i]->color.g * coEfficients[1] * val * colorAtIntersection.g;
                color.g += pointLights[i]->color.g * coEfficients[2] * pow(phong,shine) * colorAtIntersection.g;

                color.b += pointLights[i]->color.b * coEfficients[1] * val * colorAtIntersection.b;
                color.b += pointLights[i]->color.b * coEfficients[2] * pow(phong,shine) * colorAtIntersection.b;

            }
        }

         for(int i = 0; i < spotlights.size(); i++){

            Vector3D lightPosition = spotlights[i]->pointLight.position;
            Vector3D lightDirection = intersectionPoint - lightPosition;
            lightDirection.normalize_vector();

            double dot = lightDirection.dot(&spotlights[i]->direction);
            double angle = acos(dot/(lightDirection.magnitude()*spotlights[i]->direction.magnitude())) * (180.0/pi);

            if(fabs(angle)<spotlights[i]->cutoff_angle){

                Ray lightRay = Ray(lightPosition, lightDirection);
                Ray normal = getNormal(intersectionPoint,lightRay);
                
                Ray reflection = Ray(intersectionPoint, lightRay.direction - normal.direction*2*(lightRay.direction.dot(&normal.direction)));
                
                double t2 = (intersectionPoint - lightPosition).magnitude();
                if(t2 < 1e-5) continue;
                
                bool obscured = false;
                
                for(Object *obj : objects){
                    double t3 = obj->objectRayIntersection(lightRay );
                    if(t3 > 0 && t3 + 1e-5 < t2){
                        obscured = true;
                        break;
                    }
                }
                
                if(!obscured){
                    
                    double phong = max(0.0,-ray.direction.dot(&reflection.direction));
                    double val = max(0.0, -lightRay.direction.dot(&normal.direction));
                    
                    color.r += spotlights[i]->pointLight.color.r * coEfficients[1] * val * colorAtIntersection.r;
                    color.r += spotlights[i]->pointLight.color.r * coEfficients[2] * pow(phong,shine) * colorAtIntersection.r;
                    
                    color.g += spotlights[i]->pointLight.color.g * coEfficients[1] * val * colorAtIntersection.g;
                    color.g += spotlights[i]->pointLight.color.g * coEfficients[2] * pow(phong,shine) * colorAtIntersection.g;
                    
                    color.b += spotlights[i]->pointLight.color.b * coEfficients[1] * val * colorAtIntersection.b;
                    color.b += spotlights[i]->pointLight.color.b * coEfficients[2] * pow(phong,shine) * colorAtIntersection.b;
                    
                }
            }
        }

       if(level < recursionLevel){
            // if(level > 1) cout << "Recursion level " << level << endl;

            // find normal at intersectionPoint
            Ray normal = getNormal(intersectionPoint,ray);

            // find reflected ray
            Ray reflectionRay = Ray(intersectionPoint, ray.direction - normal.direction*2*(ray.direction.dot(&normal.direction)));
 
            reflectionRay.origin = reflectionRay.origin + reflectionRay.direction*1e-5;
            

             // find nearest intersection object and do recursive call

            int nearestObjectIndex = -1;
            double t = -1,tMin = 1e9;

            for(int k=0;k<(int)objects.size();k++)
            {
                t = objects[k]->intersect(reflectionRay,color, 0);
                if(t> 0 && t<tMin)
                    tMin = t , nearestObjectIndex = k;
            }

            if(nearestObjectIndex != -1)
            {
                // cout<<"Object "<<nearestObjectIndex<<" intersected"<<endl;

                Color colorTemp(0,0,0); // refelction color
                // cout<<"Before Color "<<color.r<<" "<<color.g<<" "<<color.b<<endl;
                double t = objects[nearestObjectIndex]->intersect(reflectionRay,colorTemp, level+1);

                // colorTemp will be updated while in the subsequent call
                // update color using the impact of reflection
                
                color.r += colorTemp.r * coEfficients[3];
                color.g += colorTemp.g * coEfficients[3];
                color.b += colorTemp.b * coEfficients[3];

            }
            
            
            // Vector3D reflection = lightDirection - 2*(lightDirection*normal)*normal;
            // reflection.normalize();
            // double diffuse = max(0.0, lightDirection*normal);
            // double specular = pow(max(0.0, reflection*ray.dir), shine);
            // color.r += colorAtIntersection.r * coEfficients[1] * diffuse + colorAtIntersection.r * coEfficients[2] * specular;
            // color.g += colorAtIntersection.g * coEfficients[1] * diffuse + colorAtIntersection.g * coEfficients[2] * specular;
            // color.b += colorAtIntersection.b * coEfficients[1] * diffuse + colorAtIntersection.b * coEfficients[2] * specular;
        }

        return t;
    }
*/
    virtual void print()=0;

    virtual ~Object(){}
    
};

 
class Sphere: public Object{
    
    void drawSphereUsingStacks(int stacks, int slices, double radius){
        // stack: no of divisoin along y-axis(vertical axis)
        // slice: no of division along x-axis(horizontal axis)
        for (int i = 0; i <= stacks; ++i) {
            // draw a series of quadrilateral strips to form the sphere
            glBegin(GL_QUAD_STRIP);

            for (int j = 0; j <= slices; ++j) {
                float phi = i * pi / stacks;
                float theta = j * 2 * pi / slices;

                glColor3f(color.r, color.g, color.b);

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

    public:
    Sphere(){}
    Sphere(Vector3D  center, double radius){
        this->reference_point = center;
        this->length = radius;
    }

    virtual void draw(){
        glPushMatrix();
        {
            glTranslatef(reference_point.x, reference_point.y, reference_point.z);
            drawSphereUsingStacks(50,50,length);
            // glutSolidSphere(length, 50, 50);
        
        }glPopMatrix();
    }

    virtual double objectRayIntersection(Ray &ray){
        ray.origin = ray.origin - reference_point;
        
        /* 
            From slide------ 
            Quadratic: at^2 + bt + c = 0
            a = 1 (remember, ||Rd|| = 1)
            b = 2Rd·Ro
            c = Ro·Ro - r^2

            descriminant d = b^2 - 4ac
            if d<0, no intersection return t=-1
            else if d=0 and fabs(a)<epsilon , then return (t = -c/b).
            else 
                t = -b +- sqrt(d) / 2a 

            ensures t2 is greater
            check if either t1 or t2 is positive and less than or equal to 1:
                if yes, it returns the positive t value
            if both are negative or greater than 1, it returns -1 indicating no intersection.
        */ 
        
        double a = 1; 
        double b = 2 * ray.direction.dot(&ray.origin);
        double c = ray.origin.dot(&ray.origin) - length*length;

        double t = -1;
        double d = b*b - 4*a*c;


        if(d<0){
            t=-1;
        }
        else{
            if(fabs(a)<epsilon ){
                t = -c/b;
                return t;
            }
        
            double t1 = (-b - sqrt(d))/(2*a);
            double t2 = (-b + sqrt(d))/(2*a);
            
            // make t2 larger
            if(t1>t2) swap(t1, t2);
            
            if(t1>0){
                t = t1;
            }
            else if(t2>0){
                t = t2;
            }
            else{
                t = -1;
            }
        }

        return t;
    }

    

    virtual Ray getNormal(Vector3D  point, Ray ray){
        Vector3D  normal = point - reference_point;
        normal.normalize_vector();
        return Ray(point, normal);
    }

    // take sphere input from file
    friend istream &operator >>(istream &in, Sphere &s){
        in>>s.reference_point;
        in>>s.length;
        in>>s.color;
        in>>s.coEfficients[0]>>s.coEfficients[1]>>s.coEfficients[2]>>s.coEfficients[3];
        in>>s.shine;
        return in;
    }

    virtual void print(){
        cout<<"Sphere\n";
        cout<<"reference_point: "<<reference_point<<endl;
        cout<<"length: "<<length<<endl;
        cout<<"color: "<<color<<endl;
        cout<<"coEfficients: ";
        for(int i = 0; i<4; i++){
            cout<<coEfficients[i]<<" ";
        }
        cout<<endl;
        cout<<"shine: "<<shine<<endl;
    }

    ~Sphere(){}
};

class Triangle: public Object{
    Vector3D  a, b, c;
    
    double determinant(double mat[3][3]){
        return( 
            (mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])) - 
            (mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])) +
            (mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]))
        );
    }
    
    public:
        Triangle(){}
        Triangle(Vector3D  a, Vector3D  b, Vector3D  c){
            this->a = a;
            this->b = b;
            this->c = c;
        }

        virtual Ray getNormal(Vector3D  point, Ray ray){
            Vector3D  normal = (b - a)*(c - a);
            normal.normalize_vector();
            if(ray.direction.dot(&normal) < 0){
                normal = -normal;
            }

            return Ray(point, normal);
        }

        virtual void draw(){
            glPushMatrix();{
                glBegin(GL_TRIANGLES);
                {
                    glColor3f(color.r, color.g, color.b);
                    glVertex3f(a.x, a.y, a.z);
                    glVertex3f(b.x, b.y, b.z);
                    glVertex3f(c.x, c.y, c.z);
                }
                glEnd();
            }glPopMatrix();
        }

        virtual  double objectRayIntersection(Ray &ray){
            
            // from slide
            
            double rx = ray.origin.x, ry = ray.origin.y, rz = ray.origin.z;
            double dx = ray.direction.x, dy = ray.direction.y, dz = ray.direction.z;
            
            double A[3][3] = {
                {a.x - b.x, a.x - c.x, dx},
                {a.y - b.y, a.y - c.y, dy},
                {a.z - b.z, a.z - c.z, dz}
            };
            double B[3][3]={
                {a.x - rx, a.x - c.x, dx},
                {a.y - ry, a.y - c.y, dy},
                {a.z - rz, a.z - c.z, dz}
            };
            double C[3][3]={
                {a.x - b.x, a.x - rx, dx},
                {a.y - b.y, a.y - ry, dy},
                {a.z - b.z, a.z - rz, dz}
            };
            double T[3][3]={
                {a.x - b.x, a.x - c.x, a.x - rx},
                {a.y - b.y, a.y - c.y, a.y - ry},
                {a.z - b.z, a.z - c.z, a.z - rz}
            };

            double detA = determinant(A);
            double detB = determinant(B);
            double detC = determinant(C);
            double detT = determinant(T);

            double beta = detB/detA;
            double gamma = detC/detA;
            double t = detT/detA;

            if(beta + gamma < 1 && beta > 0 && gamma > 0 && t > 0)
                return t;
            else
                return -1;

        }
    

        // for file input
        friend istream &operator >>(istream &in, Triangle &t){
            in>>t.a>>t.b>>t.c;
            in>>t.color;
            in>>t.coEfficients[0]>>t.coEfficients[1]>>t.coEfficients[2]>>t.coEfficients[3];
            in>>t.shine;
            return in;
        }

        virtual void print(){
            cout<<"Triangle\n";
            cout<<"a: "<<a<<endl;
            cout<<"b: "<<b<<endl;
            cout<<"c: "<<c<<endl;
            cout<<"color: "<<color<<endl;
            cout<<"coEfficients: ";
            for(int i = 0; i<4; i++){
                cout<<coEfficients[i]<<" ";
            }
            cout<<endl;
            cout<<"shine: "<<shine<<endl;
        }

        ~Triangle(){}
        

};
 
class Floor: public Object{
    int noOfTiles;
  
    void drawCheckBoard(){
        for (int i = 0; i < noOfTiles; i++){
			for (int j = 0; j < noOfTiles; j++){
				
                if (((i + j) % 2) == 0) 
                    glColor3f(1, 1, 1);
				else 
                    glColor3f(0, 0, 0);

				glBegin(GL_QUADS);
                {
                    glVertex3f(reference_point.x + i * length, reference_point.y + (j + 1) * length, 0); // Bottom left
                    glVertex3f(reference_point.x + (i + 1) * length, reference_point.y + (j + 1) * length, 0); // Bottom right
                    glVertex3f(reference_point.x + (i + 1) * length, reference_point.y + j * length, 0); // Top right
                    glVertex3f(reference_point.x + i * length, reference_point.y + j * length, 0); // Top left
                }
                glEnd();
			}
		}
    } 
    
    
    public:
        Floor(){}
        Floor(double floorWidth, double tileWidth){
            noOfTiles = floorWidth/tileWidth;
            this->reference_point = Vector3D (-floorWidth/2.0, -floorWidth/2.0, 0);
            this->length = tileWidth;
        }

        virtual Color getColor(Vector3D  p){
            int x = (p.x - reference_point.x)/length;
            int y = (p.y - reference_point.y)/length;
            
            if(x<0 || x>noOfTiles || y<0 || y>noOfTiles){
                return Color(0,0,0);
            }

            if((x+y)%2 == 0){
                return Color(1,1,1);
            }
            else{
                return Color(0,0,0);
            }
        }

        virtual void draw(){
            glPushMatrix();
                drawCheckBoard();
            glPopMatrix();
        }

        virtual Ray getNormal(Vector3D  point, Ray ray){
            if(ray.direction.z > 0){
                return Ray(point, Vector3D (0,0,1));
            }
            else{
                return Ray(point, Vector3D (0,0,-1));
            } 
        }

        virtual double objectRayIntersection(Ray &ray){

            Vector3D  normal = Vector3D (0,0,1);
            double dotProduct = ray.direction.dot(&normal);
            
            // check for parallel ray
            if (round(dotProduct *100) == 0) {
                return -1;
            }

            double t = -(normal.dot(&ray.origin))/dotProduct;
            Vector3D  intersection_point = ray.origin + ray.direction * t;
            if (intersection_point.x < reference_point.x || 
                intersection_point.x > reference_point.x + length * noOfTiles || 
                intersection_point.y < reference_point.y ||
                intersection_point.y > reference_point.y + length * noOfTiles)
            {
                return -1;
            }
            return t; 
        }
      
/* 
    virtual double objectRayIntersection(Ray &ray ){
        Vector3D normal = Vector3D(0, 0, 1);
        double dotP = normal.dot(& ray.direction);
        
        if (round(dotP * 100) == 0)
			return -1;

        double t = -(normal.dot(&ray.origin)) / dotP;

        Vector3D p = ray.origin + ray.direction * t;

        if(p.x <= reference_point.x || p.x >= abs(reference_point.x) && p.y <= reference_point.y && p.y >= abs(reference_point.y)){
            return -1;
        }
        
        return t;
    }
*/
   
        
        virtual void print(){
            cout<<"Floor\n";
            cout<<"reference_point: "<<reference_point<<endl;
            cout<<"length: "<<length<<endl;
            cout<<"noOfTiles: "<<noOfTiles<<endl;
        }

        ~Floor(){}  

};

/*
class Floor : public Object{
public:
    int tiles;

    Floor(){
        tiles = 1;
    }

    Floor(int floorWidth,int tileWidth){
        tiles = floorWidth / tileWidth;
        reference_point = Vector3D(-floorWidth / 2, -floorWidth / 2, 0);
        length = tileWidth;
    }

    virtual Color getColorAt(Vector3D point){

        int tileX = (point.x - reference_point.x) / length;
		int tileY = (point.y - reference_point.y) / length;

        if(tileX<0 || tileX>=tiles || tileY<0 || tileY>=tiles){
            return Color(0,0,0);
        }

		if (((tileX + tileY) % 2) == 0)
		{
			return Color(1,1,1);
		}
		else
		{
            // cout<<"Black"<<endl;
			return Color(0,0,0);
		}
    }

    virtual Ray getNormal(Vector3D point, Ray incidentRay){
        if(incidentRay.direction.z > 0) return Ray(point, Vector3D(0, 0, 1));
        else return Ray(point, Vector3D(0, 0, -1));
    }

    virtual void draw(){
        for (int i = 0; i < tiles; i++)
		{
			for (int j = 0; j < tiles; j++)
			{
				if (((i + j) % 2) == 0) glColor3f(1, 1, 1);
				else glColor3f(0, 0, 0);

				glBegin(GL_QUADS);
				{
					glVertex3f(reference_point.x + i * length, reference_point.y + j * length, 0);
					glVertex3f(reference_point.x + (i + 1) * length, reference_point.y + j * length, 0);
					glVertex3f(reference_point.x + (i + 1) * length, reference_point.y + (j + 1) * length, 0);
					glVertex3f(reference_point.x + i * length, reference_point.y + (j + 1) * length, 0);
				}
				glEnd();
			}
		}
    }

    virtual double objectRayIntersection(Ray &ray ){
        Vector3D normal = Vector3D(0, 0, 1);
        double dotP = normal.dot(& ray.direction);
        
        if (round(dotP * 100) == 0)
			return -1;

        double t = -(normal.dot(&ray.origin)) / dotP;

        Vector3D p = ray.origin + ray.direction * t;

        if(p.x <= reference_point.x || p.x >= abs(reference_point.x) && p.y <= reference_point.y && p.y >= abs(reference_point.y)){
            return -1;
        }
        
        return t;
    }

    virtual void print(){}
    ~Floor(){}
};
*/

class General: public Object{
    //Equation: F(x,y,z) = Ax2+By2+Cz2+Dxy+Exz+Fyz+Gx+Hy+Iz+J = 0
    double A, B, C, D, E, F, G, H, I, J;

    void calculate_co_efficients(Ray &ray, double &a, double &b, double &c){
        double x0 = ray.origin.x;
        double y0 = ray.origin.y;
        double z0 = ray.origin.z;

        double xd = ray.direction.x;
        double yd = ray.direction.y;
        double zd = ray.direction.z;

        a = A*xd*xd + B*yd*yd + C*zd*zd + D*xd*yd + E*xd*zd + F*yd*zd;
        b = 2*A*x0*xd + 2*B*y0*yd + 2*C*z0*zd + D*(x0*yd + y0*xd) + 
            E*(x0*zd + z0*xd) + F*(y0*zd + z0*yd) + G*xd + H*yd + I*zd;
        c = A*x0*x0 + B*y0*y0 + C*z0*z0 + D*x0*y0 + E*x0*z0 + 
            F*y0*z0 + G*x0 + H*y0 + I*z0 + J;

    }

    bool solve_quadratic(double a, double b, double c, double &t1, double &t2){
        double descriminant = b*b - 4*a*c;
        if(descriminant < 0) return false;
        
        t1 = (-b - sqrt(descriminant))/(2*a);
        t2 = (-b + sqrt(descriminant))/(2*a);
        return true;
    }

    bool check(Vector3D p){
        if(fabs(length) > epsilon){
            if(p.x < reference_point.x) return false;
            if(p.x > reference_point.x + length) return false;
        }
        if(fabs(width) > epsilon){
            if(p.y < reference_point.y) return false;
            if(p.y > reference_point.y + width) return false;
        }
        if(fabs(height) > epsilon){
            if(p.z < reference_point.z) return false;
            if(p.z > reference_point.z + height) return false;
        }
        return true;
    }

    public:
        General(){}
        General(double A, double B, double C, double D, double E, 
        double F, double G, double H, double I, double J): A(A), B(B), C(C), D(D), E(E), 
                                                           F(F), G(G), H(H), I(I), J(J){}
        
        virtual void print(){
            cout<<"General\n";
            cout<<"A: "<<A<<endl;
            cout<<"B: "<<B<<endl;
            cout<<"C: "<<C<<endl;
            cout<<"D: "<<D<<endl;
            cout<<"E: "<<E<<endl;
            cout<<"F: "<<F<<endl;
            cout<<"G: "<<G<<endl;
            cout<<"H: "<<H<<endl;
            cout<<"I: "<<I<<endl;
            cout<<"J: "<<J<<endl;
        }

        virtual Ray getNormal(Vector3D  point, Ray ray){
            double x = point.x, y = point.y, z = point.z;
            double nx = 2*A*x + D*y + E*z + G;
            double ny = 2*B*y + D*x + F*z + H;
            double nz = 2*C*z + E*x + F*y + I;
            Vector3D  normal = Vector3D (nx, ny, nz);
            normal.normalize_vector();
            return Ray(point, normal);
        }

        virtual void draw(){
            // no need to draw
            return;
        }

        virtual double objectRayIntersection(Ray &ray){
            double a, b, c;
            /*
                why a, b and c ?
                F(x,y,z) = Ax2+By2+Cz2+Dxy+Exz+Fyz+Gx+Hy+Iz+J = 0
                putting the parametric form of ray
                    F(x0+xd*t, y0+yd*t, z0+zd*t) = 0
                after calculation we get
                    a*t^2 + b*t + c = 0
            */
            
            calculate_co_efficients(ray, a, b, c);

            double t1, t2;
            
            bool isIntersecting = solve_quadratic(a, b, c, t1, t2);
            
            if(!isIntersecting) return -1;
            else{
                if(t1<0 && t2<0) return -1;

                if(t1>t2) swap(t1, t2); // t2 is always greater

                if(t1>0){
                    Vector3D  intersection_point = ray.origin + ray.direction*t1;
                    if(check(intersection_point)){
                        return t1;
                    }
                }
                if(t2 > 0){
                    Vector3D  intersection_point = ray.origin + ray.direction*t2;
                    if(check(intersection_point)){
                        return t2;
                    }
                }
                return -1;
            }

        }

        friend istream &operator >>(istream &in, General &obj){
            in>>obj.A>>obj.B>>obj.C>>obj.D>>obj.E>>obj.F>>obj.G>>obj.H>>obj.I>>obj.J;
            in>>obj.reference_point>>obj.length>>obj.width>>obj.height;
            in>>obj.color;
            in>>obj.coEfficients[0]>>obj.coEfficients[1]>>obj.coEfficients[2]>>obj.coEfficients[3];
            in>>obj.shine;
            return in;
        }

        ~General(){}    
};
 

/*
class General : public Object{
public:
    double A,B,C,D,E,F,G,H,I,J;

    General(){

    }

    virtual void print(){}

    virtual void draw(){
        return;
    }

    virtual Ray getNormal(Vector3D point, Ray incidentRay)
    {
        // Vector3D dir(2*A*point.x + D*point.y + E*point.z + G,
        //        2*B*point.y + D*point.x + F*point.z + H,
        //        2*C*point.z + E*point.x + F*point.y + I);

        // return Ray(point, dir);
        // virtual Ray getNormal(Vector3D  point, Ray ray){
            double x = point.x, y = point.y, z = point.z;
            double nx = 2*A*x + D*y + E*z + G;
            double ny = 2*B*y + D*x + F*z + H;
            double nz = 2*C*z + E*x + F*y + I;
            Vector3D  normal = Vector3D (nx, ny, nz);
            normal.normalize_vector();
            return Ray(point, normal);
        // }
    }

    bool ok(Vector3D point)
    {
        if(fabs(length) > 1e-5){
            if(point.x < reference_point.x) return false;
            if(point.x > reference_point.x + length) return false;
        }
        

        if(fabs(width) > 1e-5){
            if(point.y < reference_point.y) return false;
            if(point.y > reference_point.y + width) return false;
        }
        

        if(fabs(height) > 1e-5){
            if(point.z < reference_point.z) return false;
            if(point.z > reference_point.z + height) return false;
        }
    
        return true;
    }


    virtual double objectRayIntersection(Ray &ray){

        double X0 = ray.origin.x;
        double Y0 = ray.origin.y;
        double Z0 = ray.origin.z;

        double X1 = ray.direction.x;
        double Y1 = ray.direction.y;
        double Z1 = ray.direction.z;

        double C0 = A*X1*X1 + B*Y1*Y1 + C*Z1*Z1 + D*X1*Y1 + E*X1*Z1 + F*Y1*Z1;
        double C1 = 2*A*X0*X1 + 2*B*Y0*Y1 + 2*C*Z0*Z1 + D*(X0*Y1 + X1*Y0) + E*(X0*Z1 + X1*Z0) + F*(Y0*Z1 + Y1*Z0) + G*X1 + H*Y1 + I*Z1;
        double C2 = A*X0*X0 + B*Y0*Y0 + C*Z0*Z0 + D*X0*Y0 + E*X0*Z0 + F*Y0*Z0 + G*X0 + H*Y0 + I*Z0 + J;

        double discriminant = C1*C1 - 4*C0*C2;
        if(discriminant < 0) return -1;
        if(fabs(C0) < 1e-5) {
            return -C2/C1;
        }
        double t1 = (-C1 - sqrt(discriminant))/(2*C0);
        double t2 = (-C1 + sqrt(discriminant))/(2*C0);

        if(t1 < 0 && t2 < 0) return -1;

        // cout<<"t1 "<<t1<<" t2 "<<t2<<endl;

        if(t2<t1) swap(t1,t2);

        if(t1 > 0) {
            // cout<<"t1 "<<t1<<endl;
            Vector3D intersectionPoint = ray.origin + ray.direction*t1;
            if(ok(intersectionPoint)){
                return t1;
            }
        }
        if(t2 > 0) {
            // cout<<"t2 "<<t2<<endl;
            Vector3D intersectionPoint = ray.origin + ray.direction*t2;
            if(ok(intersectionPoint)){
                return t2;
            }
        }

        return -1;

    }
    
    // input stream
    friend istream& operator>>(istream &in, General &g)
    {
        in >> g.A >> g.B >> g.C >> g.D >> g.E >> g.F >> g.G >> g.H >> g.I >> g.J;
        in >> g.reference_point >> g.length >> g.width >> g.height;

        in >> g.color.r >> g.color.g >> g.color.b; // color
        for(int i = 0; i < 4; i++) in >> g.coEfficients[i];
        in >> g.shine;
        return in;
    }

};

*/

#endif