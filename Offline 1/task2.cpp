#include "Vector3D.h"
#include "Camera.h"
#include <GL/glut.h>

#define pi (2*acos(0.0))


double edge_len ;  // length of an edge of a triangle
double max_edge_len;    // maximum length of an edge
double max_radius;    // maximum radius of a sphere, depends on edgeLength
double current_radius;    // current sphere radius, initially set to 0 and incremented 
double increment_radius;        // measure to increment the current sphere radius

Camera* camera;


void drawTriangle(double l){
    glBegin(GL_TRIANGLES);{
        glVertex3f(l,0,0);
        glVertex3f(0,l,0);
        glVertex3f(0,0,l);
    }glEnd();
}

void drawTriangles(){
	double d = max_edge_len - edge_len;
	d /=3.0;


    for(int i = 0; i<4; i++){
        glPushMatrix();{
            glColor3f((i+1)%2, i%2, 1);
            glRotatef(i*90, 0, 1, 0);
			glTranslatef(d, d, d);
            drawTriangle(edge_len);
        }glPopMatrix(); 
    }

    for(int i = 0; i<4; i++){
        glPushMatrix();{
            glColor3f((i)%2, (i+1)%2, 1);
            glRotatef(i*90, 0, 1, 0);
            glRotatef(180, 1, 0, 1);
			glTranslatef(d,d,d);
            drawTriangle(edge_len);
        }glPopMatrix(); 
    }
}

void drawCylinderPart(double radius, double height, int segments){
	Vector3D points[segments+1];
	double phi = 70.5287794 * pi/180;		// in radians

	for(int i = 0; i<segments+1; i++){
		double theta = -phi/2 + i * phi/(segments-1);  // spread the points from -phi/2 to phi/2
		points[i].x = radius * cos(theta);
		points[i].y = radius * sin(theta);
	}

	glBegin(GL_QUADS);{
		for(int i = 0; i<segments; i++){
			glVertex3f(points[i].x, points[i].y, height/2);
			glVertex3f(points[i].x, points[i].y, -height/2);
			glVertex3f(points[i+1].x, points[i+1].y, -height/2);
			glVertex3f(points[i+1].x, points[i+1].y, height/2);
		}
	}glEnd();

}

void drawSpherePart(double radius,int slices )
{
	struct Vector3D points[slices+1][slices+1];
	//generate points
	for(int i=0;i<=slices;i++)
	{
		for(int j=0;j<=slices;j++)
		{
			// create a plane above the xy plane at point z = 1; 
			points[i][j].x= -1 +  (double)i/slices*2 ;
			points[i][j].y= -1 +  (double)j/slices*2 ;
			points[i][j].z=1;

			// because of normalization, the points are now on the surface of a sphere
			points[i][j] = points[i][j].normalize_vector();
			// now multiply the points with the radius to get the desired radius	
			points[i][j] = points[i][j]*radius;
		}
	}
	//draw quads using generated points
	for(int i=0;i<slices;i++)
	{
        // glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(int j=0;j<slices;j++)
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

void drawSphere()
{
	for(int i=0;i<4;i++){
        glPushMatrix();{
            glColor3f(0, (i+1)%2, i%2);     
            glRotatef(90*i,0,1,0);
            glTranslatef(0,0,edge_len);
            drawSpherePart(current_radius,100);
        }glPopMatrix();
    }

	for(int i=0;i<2;i++){
        glPushMatrix();{
            glColor3f(1,0,0);     
            glRotatef(90+180*i,1,0,0);
            glTranslatef(0,0,edge_len);
            drawSpherePart(current_radius,100);
        }glPopMatrix();
    }
}

void drawCylinder(){
	
	// 8 cylinder parts for 8 edges of the octahedorn

	//upper part
	for(int i = 0; i<4; i++){
		glPushMatrix();{
			glRotatef(45+90*i, 0, 1, 0);
			glTranslatef(edge_len/sqrt(2), 0, 0);
			drawCylinderPart(current_radius, edge_len*sqrt(2), 100);	
		}glPopMatrix();
	}

	// lower part
	for(int i = 0; i<4; i++){
		glPushMatrix();{
			glRotatef(90, 0, 0, 1);
			glRotatef(45+90*i, 0, 1, 0);
			glTranslatef(edge_len/sqrt(2), 0, 0);
			drawCylinderPart(current_radius, edge_len*sqrt(2), 100);	
		}glPopMatrix();
	}

	// on the xy plane
	for(int i = 0; i<4; i++){
		glPushMatrix();{	
			glRotatef(45+90*i, 0, 0, 1);
			glRotatef(90, 1, 0, 0);
			glTranslatef(edge_len/sqrt(2), 0, 0);
			drawCylinderPart(current_radius, edge_len*sqrt(2), 100);
		}glPopMatrix();
	}
}


void keyboardListener(unsigned char key, int x,int y){
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
		case '.':
			edge_len += 0.1;
			current_radius -= increment_radius;
			if(edge_len > max_edge_len){
				edge_len = max_edge_len;
				current_radius = 0.0;
			}
			break;
		case ',':
			edge_len -= 0.1;
			current_radius += increment_radius;
			if(current_radius > max_radius){
				current_radius = max_radius;
				edge_len = 0.0;
			}
			break;

		case ' ':
			cout<<"Current Camera position";camera->position->print_vector();
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

    
	drawTriangles();
	drawSphere();
	drawCylinder();

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate(){
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){

	edge_len = 1.6;
	max_edge_len = 1.6;
    
	//The formula for the radius of the inscribed sphere\
    in a regular tetrahedron is edgeLength / sqrt(3.0)
    max_radius = edge_len / sqrt(3.0);
    current_radius = 0.0;
    increment_radius = max_radius / 16;        // 16 steps to reach the maximum radius
    
    camera = new Camera(
        new Vector3D(-0.127056, -4.212026, 2.277302),
        new Vector3D(0, 0, 0),
        new Vector3D(0, 1, 0)
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