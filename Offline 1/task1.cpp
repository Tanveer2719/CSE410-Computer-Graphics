#include "Vector3D.h"
#include "Camera.h"
#include "Sphere.h"

#include <GL/glut.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
double angle;
bool simulationMode;

Camera* camera;
Sphere* sphere;

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

void drawCheckBoard(){
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

    tile_s = .1;
    int small_box_length = 51, small_box_width = 51;
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
		case ' ':
			simulationMode = !simulationMode;
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

void timerHandler(int value) {
	if (simulationMode) {
		sphere->move_forward();
	}
	glutTimerFunc(100, timerHandler, 0);	// recursive call after 100 miliseconds
}

void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	gluLookAt(camera->position->x, camera->position->y, camera->position->z,
              camera->center->x, camera->center->y, camera->center->z,
              camera->up->x, camera->up->y, camera->up->z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


    glRotatef(45,0,0,1);
	drawCheckBoard();  
 
	sphere->drawAll();

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
	simulationMode = false;

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

	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);

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
	glutTimerFunc(200, timerHandler, 0); // timerHandler is called after 200 milliseconds

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}