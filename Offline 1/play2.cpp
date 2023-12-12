#include "Vector3D.h"
#include "Camera.h"
#include "Sphere.h"

#include <GL/glut.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

Camera* camera;
Sphere* sphere;


struct point
{
	double x,y,z;
};


void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}

void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,0);
		glVertex3f( a,-a,0);
		glVertex3f(-a,-a,0);
		glVertex3f(-a, a,0);
	}glEnd();
}

void drawFloor1(){
	double tile_s = 1.5;
	int length  = 50, width = 50;
	bool white = false;

    // make length and width odd
    length += ((1+length) % 2);
    width += ((1+width) % 2);

	glPushMatrix();{	
		glTranslatef(-length * tile_s, -width * tile_s, 0);
		// glRotatef(1.5, 0, 0, 1);
		for(int i = 0; i<length; i++){
			for (int j = 0; j<width; j++){
				if(white){
					glColor3f(1, 1, 1);		// color white
					white = false;
				}
				else{
					glColor3f(0, 0, 0); 	// color black
					white = true;
				}

				glPushMatrix();
				{	glTranslatef(2*j*tile_s, 2*i*tile_s, 0);
					drawSquare(tile_s);
				}
				glPopMatrix();
			}
		}
	} glPopMatrix();

    tile_s = .4;
    int small_box_length = 21, small_box_width = 21;
    glPushMatrix();{
        glTranslatef(-small_box_length * tile_s, -small_box_width * tile_s, 0);
        for(int i = 0; i<small_box_length; i++){
            for (int j = 0; j<small_box_width;){
                
                glColor3f(1, 0, 0);     // color red

                glPushMatrix();
                {   glTranslatef(2*j*tile_s, 2*i*tile_s, 1);
                    drawSquare(tile_s);
                }
                glPopMatrix();

                if(j == small_box_length-1) break;
                if (i == 0 || i == small_box_width - 1) j++;
                else j = small_box_length - 1;
            }
        }
    } glPopMatrix();

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
		case 'i':
			sphere->move_forward();
			break;
		case 'k':
			sphere->move_backward();
			break;
		case 'j':
			sphere->rotate_arrow_anticlockwise();
			break;
		case 'l':
			sphere->rotate_arrow_clockwise();
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

void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
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

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(camera->position->x, camera->position->y, camera->position->z,
              camera->center->x, camera->center->y, camera->center->z,
              camera->up->x, camera->up->y, camera->up->z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	// drawAxes();
	// drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    // drawSS();

    glRotatef(45,0,0,1);
	drawFloor1();  
 

    // drawSquare(10);
    // sphere->drawSphere();
	// sphere->drawArrow();

	// glRotatef(-45,0,0,1);
	sphere->drawAll();

    //drawCircle(30,24);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;

    camera = new Camera(
        new Vector3D(0, -3.706431,  4.728740),
        new Vector3D(0, 0, 0),
        new Vector3D(0, 1, 0)
    );

    sphere = new Sphere(
        new Vector3D(0, 0, 0),
        new Vector3D(0, 1, 0),
        new Vector3D(0, 0, 0)
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
	glutInitWindowSize(600, 600);
	glutInitWindowPosition(400, 200);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}