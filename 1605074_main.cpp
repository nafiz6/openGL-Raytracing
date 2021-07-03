#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>

#include <windows.h>
#include <GL/glut.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include "1605074_classes.h"
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))

using namespace std;


double windowWidth = 500;
double windowHeight = 500;
double imageWidth, imageHeight;


vector <Object*> objects;
vector <Light> lights;
bitmap_image image;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double viewAngle = 80;



class Camera {
	int mag = 2;
	public:
	 Vector3D position;
	 Vector3D lookPosition;
	 Vector3D upVector;
	 Vector3D lookVector;

	 void updateCamera(){
		 lookPosition.x = position.x + lookVector.x;
		 lookPosition.y = position.y + lookVector.y;
		 lookPosition.z = position.z + lookVector.z;
	 }

	 Camera(){
		 position.x = 200;
		 position.y = 0;
		 position.z = 10;

		 lookVector.x = -200;
		 lookVector.y = -0;
		 lookVector.z = -10;

		 upVector.x = 0;
		 upVector.y = 0;
		 upVector.z = 1;
		 updateCamera();
	 }

	 void moveForward(int dir){
		 Vector3D normalisedLook = normalize(lookVector);
		 dir*=mag;

		 position.x -= normalisedLook.x*dir;
		 position.y -= normalisedLook.y*dir;
		 position.z -= normalisedLook.z*dir;

		 updateCamera();
	 }

	 void moveLR(int dir){
		 Vector3D lr = normalize(crossProduct(lookVector, upVector));
		 dir*=mag;
		 position.x += lr.x * dir;
		 position.y += lr.y * dir;
		 position.z += lr.z * dir;

		 updateCamera();
	 }

	 void moveUD(int dir){
		 dir*=mag;
		 position.x += upVector.x * dir;
		 position.y += upVector.y * dir;
		 position.z += upVector.z * dir;

		 updateCamera();
	 }

	 void rotateLR(int dir){
		 double degree = -3 * dir * pi /180.0;
		 rotate(&lookVector, upVector, degree);

		 updateCamera();
	 }

	 void rotateUD(int dir){
		 double degree = -3 * dir * pi /180.0;
		 Vector3D right = normalize(crossProduct(lookVector, upVector));
		 rotate(&lookVector, right, degree);
		 rotate(&upVector, right, degree);

		 updateCamera();
	 }

	 void tiltCW(int dir){
		 double degree = -3 * dir * pi /180.0;
		 rotate(&upVector, lookVector, degree);

		 updateCamera();
	 }
};


Camera camera;
vector<Vector3D> bullets;

Vector3D bulletCollision(){
	Vector3D res;
//	res.x = tan((-cylinder.angleY - rBall.angleY)*pi / 180.0) * (800 - rBall.radius);
	//res.y = tan((cylinder.angleX + rBall.angleX)*pi / 180.0) * (800 - rBall.radius);
	res.z = 800 - 3;

	//cout << "BULLET " << res.x << " " << res.y << endl;
	return res;
}


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


void drawCenterSphere(int type, double radius){
	int slices = 60;
	int stacks = 20;
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
		for(j=0;j<slices;j++)
		{
			if (j%2 == 0)
				glColor3f(1.0,1.0,1.0);
			else
				glColor3f(0,0,0);

			glBegin(GL_QUADS);{
			    //upper hemisphere
				if (type == 1){
					glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
					glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
					glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
					glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
				}
                //lower hemisphere
				
				else {
					glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
					glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
					glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
					glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
				}
			}glEnd();
		}
	}


}

void drawCenterConvex(double radius){
	int slices = 60;
	int stacks = 20;
	Vector3D points[100][100];
	int i,j;
	double h,r,r2;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=(2*radius - r)*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=(2*radius - r)*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
		for(j=0;j<slices;j++)
		{
			if (j%2 == 0)
				glColor3f(1.0,1.0,1.0);
			else
				glColor3f(0,0,0);

			glBegin(GL_QUADS);{
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
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



void capture(){
	cout << "Capturing\n" << flush;
    // initialize as black
    for (int i = 0; i < imageWidth; i++)
    {
        for (int j = 0; j < imageHeight; j++)
        {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }


	double planeDistance = (windowHeight / 2.0) / tan((pi * viewAngle) / (2.0*180));

	Vector3D rightVector = normalize(crossProduct(camera.lookVector, camera.upVector));
	Vector3D topLeft = vectorAdd(camera.position,
								vectorAdd( 
									vectorSub(scale(normalize(camera.lookVector), planeDistance), 
											scale(rightVector, windowWidth / 2.0)), 
									scale(normalize(camera.upVector), windowHeight / 2.0)
								)
	);

	double du = windowWidth / imageWidth;
	double dv = windowHeight / imageHeight;

	topLeft = vectorAdd(
				topLeft, 
				vectorSub(
					scale(rightVector, 0.5 * du),
					scale(camera.upVector, 0.5 * dv)
				)
			);
	
	double t, tMin;

	for (int i=1; i<imageWidth; i++){
		for (int j=1; j<imageHeight; j++){
			Vector3D curPixel = vectorAdd(
				topLeft,
				vectorSub(
					scale(rightVector, i * du),
					scale(camera.upVector, j * dv)
				)
			);
			Ray ray = Ray(camera.position, vectorSub(curPixel, camera.position));
			double *color = new double[3];
			color[0] = 1;
			color[1] = 0;
			color[2] = 0;
			float minT = -1;
			int nearest=-1;
			for (int k=0; k<objects.size(); k++){
				double *dummyColor = new double[3];
				float t = objects[k]->intersect(ray, dummyColor, 0);
				if ((nearest == -1 ||  t < minT) && t != -1){
					nearest = k;
					minT = t;
				}

				//upadate t so it stones min +ve val??
				// save nearest obj??

			}
			if (minT > 0){
				minT = objects[nearest]->intersect(ray, color, 1);
				//cout << color[0] << " " << color[1] << " " << color[2] << endl;
				image.set_pixel(i, j, 255 * color[0], 255 * color[1], 255 * color[2]);
			}
			//tmin = on->intersect(ray, color, 1);
			//update image pixel (i,j);

		}
	}
	cout << "capturing complete" << endl;
    image.save_image("out.bmp");


}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '0':
			capture();
			break;
		case '1':
			camera.rotateLR(1);
			break;

		case '2':
			camera.rotateLR(-1);
			break;

		case '3':
			camera.rotateUD(1);
			break;

		case '4':
			camera.rotateUD(-1);
			break;

		case '5':
			camera.tiltCW(1);
			break;

		case '6':
			camera.tiltCW(-1);
			break;

		default:
			break;
	}
	/*
	cout <<"key: " <<  key << endl;
	cout << "position: " << camera.position.x << " " << camera.position.y << " " << camera.position.z <<endl;
	cout << "look: " << camera.lookPosition.x << " " << camera.lookPosition.y << " " << camera.lookPosition.z <<endl;
	cout << "up: " << camera.upVector.x << " " << camera.upVector.y << " " << camera.upVector.z <<endl;
	cout << endl;
	*/
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			camera.moveForward(1);
			break;
		case GLUT_KEY_UP:		// up arrow key
			camera.moveForward(-1);
			break;

		case GLUT_KEY_RIGHT:
			camera.moveLR(1);
			break;
		case GLUT_KEY_LEFT:
			camera.moveLR(-1);
			break;

		case GLUT_KEY_PAGE_UP:
			camera.moveUD(1);
			break;
		case GLUT_KEY_PAGE_DOWN:
			camera.moveUD(-1);
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
	/*
	cout << "key: " << key << endl;
	cout << "position: " << camera.position.x << " " << camera.position.y << " " << camera.position.z <<endl;
	cout << "look: " << camera.lookPosition.x << " " << camera.lookPosition.y << " " << camera.lookPosition.z <<endl;
	cout << "up: " << camera.upVector.x << " " << camera.upVector.y << " " << camera.upVector.z <<endl;
	cout << endl;
	*/
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				//shoot
				Vector3D bullet = bulletCollision();	
				if (bullet.x <= 400 && bullet.x >=-400 && bullet.y <= 400 && bullet.y >= -400){
					bullets.push_back(bullet);
					//cout << "BULLET PUSHED" << endl;
				}
			}
			break;

		case GLUT_RIGHT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
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
//	gluLookAt(camera.position.x, camera.position.y, camera.position.z,	0,0,0,	0,1,0);
	gluLookAt(camera.position.x, camera.position.y, camera.position.z, 
		      camera.lookPosition.x, camera.lookPosition.y, camera.lookPosition.z,
			  camera.upVector.x, camera.upVector.y, camera.upVector.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	drawGrid();
	for (int i=0; i<objects.size(); i++){
		Object *sphere = objects[i];
		sphere->draw();

	}

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

	//drawSphere(30,24,20);




	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
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
	gluPerspective(viewAngle,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
	loadData();
	//cout << "CALLED" << std::flush;
	image = bitmap_image(imageWidth, imageHeight);

}

int main(int argc, char **argv){

	glutInit(&argc,argv);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(0, 0);
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
