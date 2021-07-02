#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include<math.h>
#include <GL/glut.h>

using namespace std;
#define pi (2*acos(0.0))

class Vector3D {
	public:
		double x,y,z;
		Vector3D(double i, double j, double k){
			x = i;
			y = j;
			z = k;
		}
		Vector3D(){};
};

class Object {
	public:
		Vector3D reference_point;
		Vector3D a, b, c;
		double height, width, length;
		double color[3];
		double coefficients[4]; // reflection coefficients
		int shine; // exponent term of specular component
		int A,B,C,D,E,F,G,H,I,J;
		Object(){
			cout << "Drawing obj\n";
		}
		virtual void draw(){
		}
		void setColor(double c[]){
			color[0] = c[0];
			color[1] = c[1];
			color[2] = c[2];
		}
		void setShine(int s){
			shine = s;
		}
		void setCoEfficients(double c[]){
			coefficients[0] = c[0];
			coefficients[1] = c[1];
			coefficients[2] = c[2];
			coefficients[3] = c[3];
		}

};

class Sphere: public Object {
	public:
		Sphere(){}
		Sphere(Vector3D center, double radius){
			reference_point = center;
			length = radius;
		}

		void draw(){
			Vector3D points[100][100];
			int i, j;
			double h, r;
			//generate points
			int stacks = 20;
			int slices = 24;
			for (i = 0; i <= stacks; i++)
			{
				h = length * sin(((double)i / (double)stacks) * (pi / 2));
				r = length * cos(((double)i / (double)stacks) * (pi / 2));
				for (j = 0; j <= slices; j++)
				{
					points[i][j].x = r * cos(((double)j / (double)slices) * 2 * pi);
					points[i][j].y = r * sin(((double)j / (double)slices) * 2 * pi);
					points[i][j].z = h;
				}
			}
			//draw quads using generated points
			for (i = 0; i < stacks; i++)
			{
				glPushMatrix();
				glTranslatef(reference_point.x, reference_point.y, reference_point.z);
				glColor3f(color[0], color[1], color[2]);
				for (j = 0; j < slices; j++)
				{
					glBegin(GL_QUADS);
					{
						//upper hemisphere
						glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
						glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
						glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
						glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
						//lower hemisphere
						glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
						glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
						glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
						glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
					}
					glEnd();
				}
				glPopMatrix();
			}
		}

};

class Triangle: public Object {
	public:
		Triangle(){}
		Triangle(Vector3D i, Vector3D j, Vector3D k){
			a = i;
			b = j;
			c = k;
		}

		void draw(){

		}

};

class Floor: public Object{
public:
	Floor(int floorWidth, int tileWidth){
		reference_point = Vector3D(-floorWidth/2, -floorWidth/2, 0);
		length = tileWidth;

	}

	void draw(){

	}

};

class Light {
public:
	Vector3D light_pos;
	double color[3];
};


// declaration
extern vector <Object*> objects;
extern vector <Light> lights;



void loadData(){
    ifstream sceneInput("scene.txt"); 
	string input;
	getline(sceneInput, input);

	int recursionLevel;
	int resolution;
	int objectCount;

	stringstream lineRecursion(input);
	lineRecursion >> recursionLevel;

	getline(sceneInput, input);
	stringstream lineRes(input);
	lineRes >> resolution;


    getline(sceneInput, input);
    getline(sceneInput, input);
	stringstream lineObjCnt(input);
	lineObjCnt >> objectCount;

	for (int i=0; i<objectCount; i++){
		string objectName;
		getline(sceneInput, objectName);

		if (objectName == "sphere"){
			Object *sphere;
			sphere = new Sphere();


			getline(sceneInput, input);
			stringstream lineCenter(input);
			lineCenter >> sphere->reference_point.x;
			lineCenter >> sphere->reference_point.y;
			lineCenter >> sphere->reference_point.z;

			getline(sceneInput, input);
			stringstream sphereRad(input);
			sphereRad >> sphere->length;

			getline(sceneInput, input);
			stringstream lineColor(input);
			lineColor >> sphere->color[0];
			lineColor >> sphere->color[1];
			lineColor >> sphere->color[2];

			getline(sceneInput, input);
			stringstream lineCoefficients(input);
			lineCoefficients >> sphere->coefficients[0];
			lineCoefficients >> sphere->coefficients[1];
			lineCoefficients >> sphere->coefficients[2];
			lineCoefficients >> sphere->coefficients[3];

			getline(sceneInput, input);
			stringstream lineShiny(input);
			lineShiny >> sphere->shine;

			objects.push_back(sphere);

		}
		else if (objectName == "triangle"){
			Object *triangle;
			triangle = new Triangle();


			getline(sceneInput, input);
			stringstream lineA(input);
			lineA >> triangle->a.x;
			lineA >> triangle->a.y;
			lineA >> triangle->a.z;

			getline(sceneInput, input);
			stringstream lineB(input);
			lineB >> triangle->b.x;
			lineB >> triangle->b.y;
			lineB >> triangle->b.z;

			getline(sceneInput, input);
			stringstream lineC(input);
			lineC >> triangle->c.x;
			lineC >> triangle->c.y;
			lineC >> triangle->c.z;


			getline(sceneInput, input);
			stringstream lineColor(input);
			lineColor >> triangle->color[0];
			lineColor >> triangle->color[1];
			lineColor >> triangle->color[2];

			getline(sceneInput, input);
			stringstream lineCoefficients(input);
			lineCoefficients >> triangle->coefficients[0];
			lineCoefficients >> triangle->coefficients[1];
			lineCoefficients >> triangle->coefficients[2];
			lineCoefficients >> triangle->coefficients[3];

			getline(sceneInput, input);
			stringstream lineShiny(input);
			lineShiny >> triangle->shine;

			objects.push_back(triangle);

		}
		else if (objectName=="general"){
			Object *object;
			object = new Object();
			getline(sceneInput, input);
			stringstream lineEqn(input);
			lineEqn >> object->A;
			lineEqn >> object->B;
			lineEqn >> object->C;
			lineEqn >> object->D;
			lineEqn >> object->E;
			lineEqn >> object->F;
			lineEqn >> object->G;
			lineEqn >> object->H;
			lineEqn >> object->I;
			lineEqn >> object->J;


			getline(sceneInput, input);
			stringstream lineRef(input);
			lineRef >> object->reference_point.x;
			lineRef >> object->reference_point.y;
			lineRef >> object->reference_point.z;

			getline(sceneInput, input);
			stringstream lineColor(input);
			lineColor >> object->color[0];
			lineColor >> object->color[1];
			lineColor >> object->color[2];

			getline(sceneInput, input);
			stringstream lineCoefficients(input);
			lineCoefficients >> object->coefficients[0];
			lineCoefficients >> object->coefficients[1];
			lineCoefficients >> object->coefficients[2];
			lineCoefficients >> object->coefficients[3];

			getline(sceneInput, input);
			stringstream lineShiny(input);
			lineShiny >> object->shine;

			objects.push_back(object);


		}
		else{
			cout << "\ninvalid object name\n";
		}
		getline(sceneInput, input);
	}
    getline(sceneInput, input);
	stringstream lineLightCnt(input);
	int lightCount;
	lineLightCnt >> lightCount;

	for (int i=0; i<lightCount; i++){
		Light light;

		getline(sceneInput, input);
		stringstream lineLightPos(input);

		lineLightPos >> light.light_pos.x;
		lineLightPos >> light.light_pos.y;
		lineLightPos >> light.light_pos.z;

		getline(sceneInput, input);
		stringstream lineLightCol(input);
		lineLightCol >> light.color[0];
		lineLightCol >> light.color[1];
		lineLightCol >> light.color[2];

		lights.push_back(light);

		
	}

	Object *floor;
	floor = new Floor(1000, 20);
	objects.push_back(floor);

}
