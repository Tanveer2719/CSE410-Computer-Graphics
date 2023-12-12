#include "Vector3D.h"
#include "Camera.h"
#include <GL/glut.h>

#define pi (2*acos(0.0))

int drawgrid;
int drawaxes;

double triangle_edge_len ;  // length of an edge of a triangle
double max_sphere_radius;    // maximum radius of a sphere, depends on edgeLength
double current_sphere_radius;    // current sphere radius, initially set to 0 and incremented 
double increment_radius;        // measure to increment the current sphere radius

Camera* camera;


void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glColor3f(1,0,0);   // x axis red colored
            glVertex3f( 100,0,0);
			glVertex3f( 0,0,0);

			glColor3f(0,1,0);   // y axis green colored
            glVertex3f(0,0,0);
			glVertex3f(0, 100,0);

			glColor3f(0,0,1);   // z axis blue colored
            glVertex3f(0,0, 100);
			glVertex3f(0,0,0);
		}glEnd();
	}
}

void drawGrid()
{
	int i;
    glColor3f(0.6, 0.6, 0.6);	//grey
    glBegin(GL_LINES);{
        for(i=-18;i<=18;i++){

            if(i==0)
                continue;	//SKIP the MAIN axes

            //lines parallel to Y-axis
            glVertex3f(i*0.1, -1.9, 0);
            glVertex3f(i*0.1,  1.9, 0);

            //lines parallel to X-axis
            glVertex3f(-1.9, i*0.1, 0);
            glVertex3f( 1.9, i*0.1, 0);
        }
    }glEnd();

}

void drawTriangle(){
    glBegin(GL_TRIANGLES);{
        glVertex3f(1,0,0);
        glVertex3f(0,1,0);
        glVertex3f(0,0,1);
    }glEnd();
}

void drawTriangles(){
    for(int i = 0; i<4; i++){
        glPushMatrix();{
            glColor3f((i+1)%2, i%2, 1);
            glRotatef(i*90, 0, 1, 0);
            drawTriangle();
        }glPopMatrix(); 
    }

    for(int i = 0; i<4; i++){
        glPushMatrix();{
            glColor3f((i)%2, (i+1)%2, 1);
            glRotatef(i*90, 0, 1, 0);
            glRotatef(180, 1, 0, 1);
            drawTriangle();
        }glPopMatrix(); 
    }
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,2);
		glVertex3f( a,-a,2);
		glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);
	}glEnd();
}

void drawCircle(double radius,int segments)
{
    int i;
    Vector3D points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
     double shade;
    Vector3D points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawSphere(double radius,int slices,int stacks)
{
	Vector3D points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}

void keyboardListener(unsigned char key, int x,int y){
	// cout<<"Current Camera position";camera->position->print_vector();
	switch(key){

		case '1':
			camera->look_left();
			break;
		case '2':
			camera->look_right();	
			break;
		case '3':
			camera->look_up();	
			break;
		case '4':
			camera->look_down();	
			break;
		case '5':
			camera->tilt_anticlokwise();
			break;
		case '6':
			camera->tilt_clockwise();
			break;
		case 'w':
			camera->up_reference();
			break;
		case 's':
			camera->down_reference();
			break;
		
		default:
			break;
	}
}

void specialKeyListener(int key, int x,int y){
	// cout<<"Current Camera position";camera->position->print_vector();
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			camera->move_backward();
			break;
		case GLUT_KEY_UP:		// up arrow key
			camera->move_forward(); 
			break;
		case GLUT_KEY_RIGHT:
			camera->move_right();
			break;
		case GLUT_KEY_LEFT:
			camera->move_left();	
			break;

		case GLUT_KEY_PAGE_UP:
			camera->move_up();
			break;
		case GLUT_KEY_PAGE_DOWN:
			camera->move_down();
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}

void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();
	
    gluLookAt(camera->position->x, camera->position->y, camera->position->z,
              camera->center->x, camera->center->y, camera->center->z,
              camera->up->x, camera->up->y, camera->up->z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/

	drawAxes();
	drawGrid();
    drawTriangles();

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=1;
	drawaxes=1;
    triangle_edge_len =  1.6;

    //The formula for the radius of the inscribed sphere\
    in a regular tetrahedron is edgeLength / sqrt(3.0)
    max_sphere_radius = triangle_edge_len / sqrt(3.0);
    current_sphere_radius = 0.0;
    increment_radius = max_sphere_radius / 16;        // 16 steps to reach the maximum radius
    
    camera = new Camera(
        new Vector3D(0, -15.907796, 9.451017),
        new Vector3D(0, 0, 0),
        new Vector3D(0, 1, 0)
    );
    

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	
    //load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	
    
    //give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(50,50);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}