#include"1905025_classes.h"
#include <GL/glut.h>


#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
bool simulationMode;
bool drawaxes = false;
bool drawgrid = false;
int no_of_image = 1;


vector<Object*> objects;
vector<PointLight*> pointLights;
vector<Spotlight*> spotlights;
int recursionLevel, image_dimension, imageHeight, imageWidth; 
int windowHeight = 500, windowWidth = 500;

int view_angle=80;

bitmap_image image;

Camera* camera;

// draw axes
void drawAxes()
{
	if(drawaxes){
		glColor3f(1.0, 1.0, 0.0);
		glBegin(GL_LINES);
		{
			glVertex3f( 500,0,0);
			glVertex3f(-500,0,0);

			glVertex3f(0,-500,0);
			glVertex3f(0, 500,0);

			glVertex3f(0,0, 500);
			glVertex3f(0,0,-500);
		}
		glEnd();
	}
}

// draw grid
void drawGrid()
{
	int i;
	if(drawgrid)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);
		{
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
		}
		glEnd();
	}
}

void loadData()
{
	ifstream input("scene.txt");
	input >> recursionLevel >> imageHeight;

	imageWidth = imageHeight;

	int no_of_objects;
	input >> no_of_objects;

	while(no_of_objects--)
	{
		string s;
		input >> s;

		Object *obj;

		if(s == "sphere"){
			obj = new Sphere();
			input >> *((Sphere *)obj);
		}
		else if(s == "triangle"){
			obj = new Triangle();
			input >> *((Triangle *)obj);
		}
		else if(s == "general"){
			// obj = new Plane();
			obj = new General();
			input >> *((General *)obj);
		}
		else{
			cout<<s<<" is not a valid object type"<<endl;
		}
		objects.push_back(obj);
	}

	int pointLightCount;
	input >> pointLightCount;

	while(pointLightCount--){
		PointLight *light = new PointLight();
		input >> *light;
		pointLights.push_back(light);
	}

	int spotlightCount;
	input >> spotlightCount;

	while(spotlightCount--){
		Spotlight *spotlight = new Spotlight();
		input >> *spotlight;
		spotlights.push_back(spotlight);
	}
    

	Object *floor = new Floor(1000, 20);
	floor->setColor(Color(0.5, 0.5, 0.5));
	floor->setcoEfficients(.4, .4, .2, .2);
	objects.push_back(floor);
	
}

void capture(){

	cout<<"inside capture\n";

	// set image background color to black
	for(int x=0; x<imageWidth; x++){
		for(int y=0; y<imageHeight; y++){
			image.set_pixel(x, y, 0, 0, 0);
		}
	}

	// plane distance is the distance from the camera to the view plane
	double plane_distance = (windowHeight/2.0) / tan(view_angle*pi/360.0);

	Vector3D topleft = *camera->position + (*camera->look*plane_distance) 
						+ (*camera->up*(windowHeight/2.0)) 
						- (*camera->right*(windowWidth/2.0)) ;

	double du = (windowWidth*1.0)/imageWidth;
	double dv = (windowHeight*1.0)/imageHeight;

	// Choose middle of the grid cell
	topleft = topleft + (*camera->right*0.5*du) - (*camera->up*0.5*dv);

	// cout<<"topleft: ";topleft.print_vector();


	
	double t, tMin;
	for(int i = 0; i<imageWidth; i++){
		for(int j = 0; j<imageHeight; j++){
			
			// current pixel
			Vector3D current = topleft + (*camera->right*i*du) - (*camera->up*j*dv);
			// cout<<"Current: ";current.print_vector();
			// break;

			
			// cast ray from camera to current pixel
			Ray ray(*camera->position, (current - *camera->position));

			Color color(0,0,0);
			tMin = -1;
			int nearest=-1 ;
			for(int k = 0; k<(int)objects.size(); k++){
				t = objects[k]->intersect(ray, color, 0);
				if(t>0 && (t<tMin || nearest == -1)){
					tMin = t;
					nearest = k;
				}
			}

			if(nearest != -1){
				color = Color(0, 0, 0);
				double t = objects[nearest]->intersect(ray, color, 1); 

				if(color.r > 1) color.r = 1;
				if(color.g > 1) color.g = 1;
				if(color.b > 1) color.b = 1;

				if(color.r<0) color.r = 0;
				if(color.g<0) color.g = 0;
				if(color.b<0) color.b = 0;

				image.set_pixel(i, j, color.r*255, color.g*255, color.b*255);
			
			}
		}
	}

	image.save_image("Output_1"+to_string(no_of_image++)+".bmp");
	cout<<"returning from capture\n";
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
		case ' ':
			capture();
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

	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	gluLookAt(camera->position->x, camera->position->y, camera->position->z,
              camera->center->x, camera->center->y, camera->center->z,
              camera->up->x, camera->up->y, camera->up->z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);
	
	drawAxes();
	// drawGrid();
	
	// draw all objects
	for(auto i: objects){
		i->draw();
	}
	for(auto i: pointLights){
		i->draw();
	}
	for (auto i: spotlights){
		i->draw();
	}

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
	simulationMode = false;
	drawaxes = true;
	drawgrid = true;

    //params: Camera(Vector3D* position, Vector3D* look, Vector3D* right)
	camera = new Camera(
        new Vector3D(210, 120,  100),
        new Vector3D(-.707, -.707, 0),
        new Vector3D(-.707, .707, 0)
    );


	loadData();
	image=bitmap_image(imageWidth, imageHeight);



	//clear the screen
	glClearColor(0,0,0,0);

	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);

}

void freeMemory()
{
	pointLights.clear();
    spotlights.clear();


    for (auto ob : objects)
    {
        delete ob;
    }
   
    objects.clear();

	cout<<"all memory freed\n";

}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(400, 200);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	atexit(freeMemory);
	glutMainLoop();		//The main loop of OpenGL

	return 0;
}