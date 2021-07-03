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

double power(double base, int exp){
	double result = 1;
	for (int i=0; i<exp; i++){
		result *= base;
	}
	return result;
}

class Light {
public:
	Vector3D light_pos;
	double color[3];

	void draw(){
		double size = 10;
		glColor3f(color[0], color[1], color[2]);

		glBegin(GL_QUADS);
		{
			glVertex3f(light_pos.x, light_pos.y, light_pos.z);
			glVertex3f(light_pos.x, light_pos.y + size, light_pos.z);
			glVertex3f(light_pos.x, light_pos.y + size, light_pos.z + size);
			glVertex3f(light_pos.x, light_pos.y, light_pos.z + size);
		}
		glEnd();
	}
};

Vector3D crossProduct(Vector3D a, Vector3D b){
	Vector3D res;
	res.x = a.y*b.z - a.z*b.y;
	res.y = -(a.x*b.z - a.z*b.x);
	res.z = a.x*b.y - a.y*b.x;
	return res;
}

float vectorDot(Vector3D a, Vector3D b){
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

Vector3D vectorSub(Vector3D a, Vector3D b){
	Vector3D res;
	res.x = a.x - b.x;
	res.y = a.y - b.y;
	res.z = a.z - b.z;
	return res;
}

Vector3D vectorAdd(Vector3D a, Vector3D b){
	Vector3D res;
	res.x = a.x + b.x;
	res.y = a.y + b.y;
	res.z = a.z + b.z;
	return res;
}

double magnitude(Vector3D a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

Vector3D normalize(Vector3D a){
	Vector3D res;
	double mag = magnitude(a);
	res.x = a.x / mag;
	res.y = a.y / mag;
	res.z = a.z / mag;
	return res;
}

Vector3D scale(Vector3D vec, double val){
	Vector3D res;
	res.x = vec.x * val;
	res.y = vec.y * val;
	res.z = vec.z * val;
	return res;
}

void rotate(Vector3D* target, Vector3D normal, double deg){
	Vector3D right = normalize(crossProduct(*target, normal));
	right = scale(right, magnitude(*target));

	target->x = target->x * cos(deg) + right.x * sin(deg);
	target->y = target->y * cos(deg) + right.y * sin(deg);
	target->z = target->z * cos(deg) + right.z * sin(deg);
}

void print(Vector3D a, string name){
	cout << name << ": " << a.x  << " " << a.y << " " << a.z << " ";
}

double determinant(double mat[3][3]){
	double a = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]);
	double b = mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]);
	double c = mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

	return a - b + c;
}


class Ray {
public:
	Vector3D start;
	Vector3D dir;

	Ray(){};
	Ray(Vector3D s, Vector3D d){
		start.x = s.x;
		start.y = s.y;
		start.z = s.z;

		dir = normalize(d);
	}

};

// Declaration
extern vector <Light> lights;
int recursionLevel;

class Object {
	public:
		Vector3D a, b, c;
		Vector3D reference_point;
		double height, width, length;
		double color[3];
		double coefficients[4]; // reflection coefficients
		int shine; // exponent term of specular component
		double A,B,C,D,E,F,G,H,I,J;

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

		virtual float intersect(Ray ray, double *color, int level);

};

// declaration
extern double imageWidth, imageHeight;
extern vector <Object*> objects;

float Object::intersect(Ray ray, double* color, int level){
	double t_2_coeff = A * (ray.dir.x * ray.dir.x) +
					   B * (ray.dir.y * ray.dir.y) +
					   C * (ray.dir.z * ray.dir.z) +
					   D * (ray.dir.x * ray.dir.y) +
					   E * (ray.dir.x * ray.dir.z) +
					   F * (ray.dir.y * ray.dir.z);

	double t_coeff = 2 * A * (ray.start.x * ray.dir.x) +
					 2 * B * (ray.start.y * ray.dir.y) +
					 2 * C * (ray.start.z * ray.dir.z) +
					 1 * D * (ray.start.x * ray.dir.y) +
					 1 * D * (ray.start.y * ray.dir.x) +
					 1 * E * (ray.start.x * ray.dir.z) +
					 1 * E * (ray.start.z * ray.dir.x) +
					 1 * F * (ray.start.y * ray.dir.z) +
					 1 * F * (ray.start.z * ray.dir.y) +
					 1 * G * (ray.dir.x) +
					 1 * H * (ray.dir.y) +
					 1 * I * (ray.dir.z);

	double constant = A * (ray.start.x * ray.start.x) +
					  B * (ray.start.y * ray.start.y) +
					  C * (ray.start.z * ray.start.z) +
					  D * (ray.start.x * ray.start.y) +
					  E * (ray.start.x * ray.start.z) +
					  F * (ray.start.y * ray.start.z) +
					  G * (ray.start.x) +
					  H * (ray.start.y) +
					  I * (ray.start.z) +
					  J;

	double descriminant = (t_coeff * t_coeff) - 4 * t_2_coeff * constant;
	if (descriminant < 0)
		return -1;

	double t1 = (-t_coeff + sqrt(descriminant)) / (2 * t_2_coeff);
	double t2 = (-t_coeff - sqrt(descriminant)) / (2 * t_2_coeff);

	Vector3D intersect1 = vectorAdd(ray.start, scale(ray.dir, t1));
	Vector3D intersect2 = vectorAdd(ray.start, scale(ray.dir, t2));

	bool t1Valid = t1 >= 0;
	bool t2Valid = t2 >= 0;
	if (round(length * 100) != 0 && abs(intersect1.x - reference_point.x) > length)
	{
		t1Valid = false;
	}
	else if (round(width * 100) != 0 && abs(intersect1.y - reference_point.y) > width)
	{
		t1Valid = false;
	}
	else if (round(height * 100) != 0 && abs(intersect1.z - reference_point.z) > height)
	{
		t1Valid = false;
	}

	if (round(length * 100) != 0 && abs(intersect2.x - reference_point.x) > length)
	{
		t2Valid = false;
	}
	else if (round(100 * width) != 0 && abs(intersect2.y - reference_point.y) > width)
	{
		t2Valid = false;
	}
	else if (round(100 * height) != 0 && abs(intersect2.z - reference_point.z) > height)
	{
		t2Valid = false;
	}

	double t = -1;
	Vector3D intersectPoint;
	if (!t1Valid && !t2Valid)
		return t;
	if (t1Valid && t2Valid)
	{
		if (t1 < t2)
		{
			intersectPoint = intersect1;
		}
		if (t1 >= t2)
		{
			intersectPoint = intersect2;
		}
		t = min(t1, t2);
	}
	else if (t1Valid)
	{
		t = t1;
		intersectPoint = intersect1;
	}
	else if (t2Valid)
	{
		t = t2;
		intersectPoint = intersect2;
	}

	if (level == 0)
		return t;

	Vector3D normal(
		2 * A * intersectPoint.x + D * intersectPoint.y + E * intersectPoint.z + G,
		2 * B * intersectPoint.y + D * intersectPoint.x + F * intersectPoint.z + H,
		2 * C * intersectPoint.z + E * intersectPoint.x + F * intersectPoint.y + I);

	normal = normalize(normal);

	double intensity[3] = {coefficients[0], coefficients[0], coefficients[0]};

	double totalSpecularIntensity[3] = {0, 0, 0};
	for (int i = 0; i < lights.size(); i++)
	{
		// Diffuse
		Vector3D L = vectorSub(intersectPoint, lights[i].light_pos);
		L = normalize(L);
		double diffuseIntensity = -1 * coefficients[1] * vectorDot(L, normal);

		// Specular
		Vector3D H = scale(vectorAdd(L, ray.dir), 0.5);
		H = normalize(H);
		double specDot = vectorDot(normal, H);
		int sign = 1;
		if (specDot < 0)
			sign = -1;
		specDot *= sign;
		double specularIntensity = -coefficients[2] * sign * power(specDot, shine);

		//cout << specularIntensity << endl;

		for (int j = 0; j < 3; j++)
		{
			intensity[j] += max(0.0, diffuseIntensity) * lights[i].color[j];
			totalSpecularIntensity[j] += max(0.0, specularIntensity) * lights[i].color[j];
		}
	}

	for (int i = 0; i < 3; i++)
	{
		color[i] = min(1.0, this->color[i] * intensity[i] + totalSpecularIntensity[i]);
		color[i] = max(color[i], 0.0);
	}

	Vector3D reflectionDir = vectorSub(ray.dir, scale(normal, 2 * vectorDot(ray.dir, normal)));
	reflectionDir = normalize(reflectionDir);
	Ray reflectRay(vectorAdd(intersectPoint, scale(reflectionDir, 0.2)), reflectionDir);
	//print(reflectionDir, "\n\nreflect dir");
	//print(intersectionPoint, "\n\nintersect point");
	//print(reference_point, "\n\nreference point");

	float minT = -1;
	int nearest = -1;
	for (int k = 0; k < objects.size(); k++)
	{
		if (this == objects[k])
			continue;
		double *dummyColor = new double[3];
		float newT = objects[k]->intersect(reflectRay, dummyColor, 0);
		if ((nearest == -1 || newT < minT) && newT > 0)
		{
			nearest = k;
			minT = newT;
		}
	}
	if (minT > 0)
	{
		double *nextColor = new double[3];
		minT = objects[nearest]->intersect(reflectRay, nextColor, level + 1);
		for (int i = 0; i < 3; i++)
		{
			//cout << nextColor[i] << " nextcolor \n";
			color[i] += nextColor[i] * coefficients[3];
			color[i] = min(color[i], 1.0);
			color[i] = max(color[i], 0.0);
		}
	}
	return t;
		}

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

		float intersect(Ray ray, double *retColor, int level){
			Vector3D newRayOrigin = vectorSub(ray.start, reference_point);
			bool outside = true;
			double dist = vectorDot(newRayOrigin, newRayOrigin);

			if (dist < length * length)
			{
				outside = false;
			}
			double tp = -vectorDot(newRayOrigin, ray.dir);

			if (outside && tp < 0)
				return -1;

			double closestDist = dist - tp * tp;
			if (closestDist > length * length)
				return -1;

			double tBar_2 = length * length - closestDist;

			double t = -1;
			if (outside)
			{
				t = tp - sqrt(tBar_2);
			}
			else
			{
				t = tp + sqrt(tBar_2);
			}

			if (level == 0)
				return t;

			Vector3D normal = vectorAdd(
				newRayOrigin,
				scale(ray.dir, t));
			normal = normalize(normal);

			Vector3D intersectionPoint = vectorAdd(
					ray.start,
					scale(ray.dir, t));

			// Diffuse
			double intensity[3] = {coefficients[0], coefficients[0], coefficients[0]};

			double totalSpecularIntensity[3] = {0, 0, 0};
			for (int i = 0; i < lights.size(); i++)
			{
				// Diffuse
				Vector3D L = vectorSub(intersectionPoint, lights[i].light_pos);
				L = normalize(L);
				double diffuseIntensity = -1 * coefficients[1] * vectorDot(L, normal);

				// Specular
				Vector3D H = scale(vectorAdd(L, ray.dir), 0.5);
				H = normalize(H);
				double specDot = vectorDot(normal, H);
				int sign = 1;
				if (specDot < 0)
					sign = -1;
				specDot *= sign;
				double specularIntensity = -coefficients[2] * sign * power(specDot, shine);

				//cout << specularIntensity << endl;

				for (int j = 0; j < 3; j++)
				{
					intensity[j] += max(0.0, diffuseIntensity) * lights[i].color[j];
					totalSpecularIntensity[j] += max(0.0, specularIntensity) * lights[i].color[j];
				}
			}

			for (int i = 0; i < 3; i++)
			{
				retColor[i] = min(1.0, this->color[i] * intensity[i] + totalSpecularIntensity[i]);
				retColor[i] = max(retColor[i], 0.0);
			}

			if (level >= recursionLevel)
				return t;

			Vector3D reflectionDir = vectorSub(ray.dir, scale(normal, 2 * vectorDot(ray.dir, normal)));
			reflectionDir = normalize(reflectionDir);
			Ray reflectRay(vectorAdd(intersectionPoint, scale(reflectionDir, 0.2)), reflectionDir);
			//print(reflectionDir, "\n\nreflect dir");
			//print(intersectionPoint, "\n\nintersect point");
			//print(reference_point, "\n\nreference point");

			float minT = -1;
			int nearest=-1;
			for (int k=0; k<objects.size(); k++){
				if (this == objects[k])continue;
				double *dummyColor = new double[3];
				float newT = objects[k]->intersect(reflectRay, dummyColor, 0);
				if ((nearest == -1 ||  newT < minT) && newT > 0){
					nearest = k;
					minT = newT;
				}

			}
			if (minT > 0){
				double *nextColor = new double[3];
				minT = objects[nearest]->intersect(reflectRay, nextColor, level+1);
					for (int i = 0; i < 3; i++)
					{
						//cout << nextColor[i] << " nextcolor \n";
						retColor[i] += nextColor[i] * coefficients[3];
						retColor[i] = min(retColor[i], 1.0);
						retColor[i] = max(retColor[i], 0.0);
					}
			}

			return t;


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

			glColor3f(color[0], color[1], color[2]);
			glBegin(GL_TRIANGLES);
			{
				glVertex3f(a.x, a.y, a.z);
				glVertex3f(b.x, b.y, b.z);
				glVertex3f(c.x, c.y, c.z);
			}
			glEnd();
		}

		float intersect(Ray ray, double *color, int level){
			double betaMat[3][3] = {
				{a.x - ray.start.x, a.x - c.x, ray.dir.x},
				{a.y - ray.start.y, a.y - c.y, ray.dir.y},
				{a.z - ray.start.z, a.z - c.z, ray.dir.z}
			};
			double gammaMat[3][3] = {
				{a.x - b.x, a.x - ray.start.x, ray.dir.x},
				{a.y - b.y, a.y - ray.start.y, ray.dir.y},
				{a.z - b.z, a.z - ray.start.z, ray.dir.z}
			};
			double tMat[3][3] = {
				{a.x - b.x, a.x - c.x, a.x - ray.start.x},
				{a.y - b.y, a.y - c.y, a.y - ray.start.y},
				{a.z - b.z, a.z - c.z, a.z - ray.start.z}
			};
			double AMat[3][3] {
				{a.x - b.x, a.x - c.x, ray.dir.x},
				{a.y - b.y, a.y - c.y, ray.dir.y},
				{a.z - b.z, a.z - c.z, ray.dir.z}
			};


			double Adet = determinant(AMat);

			double beta = determinant(betaMat) / Adet;
			double gamma = determinant(gammaMat) / Adet;
			double t = determinant(tMat) / Adet;

			//cout << Adet << beta << gamma << t << endl;
			Vector3D normal = crossProduct(
				vectorSub(b, a),
				vectorSub(c, a)
			);
			normal = normalize(normal);

			if (beta + gamma < 1 && beta > 0 && gamma > 0 && t > 0){
				if (level == 0)
					return t;
				Vector3D point = vectorAdd(
					ray.start,
					scale(ray.dir, t));
				double intensity[3] = {coefficients[0], coefficients[0], coefficients[0]};

				double totalSpecularIntensity[3] = {0, 0, 0};
				for (int i = 0; i < lights.size(); i++)
				{
					// Diffuse
					Vector3D L = vectorSub(point, lights[i].light_pos);
					L = normalize(L);
					double diffuseIntensity = -1 * coefficients[1] * vectorDot(L, normal);

					// Specular
					Vector3D H = scale(vectorAdd(L, ray.dir), 0.5);
					H = normalize(H);
					double specDot = vectorDot(normal, H);
					int sign = 1;
					if (specDot < 0)
						sign = -1;
					specDot *= sign;
					double specularIntensity = -coefficients[2] * sign * power(specDot, shine);

					for (int j = 0; j < 3; j++)
					{
						intensity[j] += max(0.0, diffuseIntensity) * lights[i].color[j];
						totalSpecularIntensity[j] += max(0.0, specularIntensity) * lights[i].color[j];
					}
				}

				for (int i = 0; i < 3; i++)
				{
					color[i] = min(1.0, this->color[i] * intensity[i] + totalSpecularIntensity[i]);
					color[i] = max(color[i], 0.0);
				}

				if (level >= recursionLevel)
					return t;

				Vector3D reflectionDir = vectorSub(ray.dir, scale(normal, 2 * vectorDot(ray.dir, normal)));
				reflectionDir = normalize(reflectionDir);
				Ray reflectRay(vectorAdd(point, scale(reflectionDir, 0.2)), reflectionDir);
				//print(reflectionDir, "\n\nreflect dir");
				//print(intersectionPoint, "\n\nintersect point");
				//print(reference_point, "\n\nreference point");

				float minT = -1;
				int nearest = -1;
				for (int k = 0; k < objects.size(); k++)
				{
					if (this == objects[k])
						continue;
					double *dummyColor = new double[3];
					float newT = objects[k]->intersect(reflectRay, dummyColor, 0);
					if ((nearest == -1 || newT < minT) && newT > 0)
					{
						nearest = k;
						minT = newT;
					}
				}
				if (minT > 0)
				{
					double *nextColor = new double[3];
					minT = objects[nearest]->intersect(reflectRay, nextColor, level + 1);
					for (int i = 0; i < 3; i++)
					{
						//cout << nextColor[i] << " nextcolor \n";
						color[i] += nextColor[i] * coefficients[3];
						color[i] = min(color[i], 1.0);
						color[i] = max(color[i], 0.0);
					}
				}

				return t;
			}
			return -1;



		}

};

class Floor: public Object{
public:
	int tileCount;
	Floor(int floorWidth, int tileWidth)
	{
		reference_point = Vector3D(-floorWidth / 2, -floorWidth / 2, 0);
		length = tileWidth;
		tileCount = floorWidth / tileWidth;
	}

	void draw()
	{
		for (int i = 0; i < tileCount; i++)
		{
			for (int j = 0; j < tileCount; j++)
			{
				if (((i + j) % 2) == 0)
				{
					glColor3f(1, 1, 1);
				}
				else
				{
					glColor3f(0, 0, 0);
				}
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
	float intersect(Ray ray, double *color, int level)
	{
		Vector3D normal = Vector3D(0, 0, 1);
		float div = vectorDot(normal, ray.dir);
		if (round(div * 100) == 0)
			return -1;
		float t = -(vectorDot(normal, ray.start)) / div;
		//cout << "T" <<  t << endl;

		Vector3D intersectionPoint = vectorAdd(ray.start, scale(ray.dir, t));

		//print(intersectionPoint, "intersectionPoint\n\n");
		//print(ray.dir, "rayDir\n\n");
		if ((intersectionPoint.x) >= abs(reference_point.x) || intersectionPoint.x <= reference_point.x || (intersectionPoint.y) >= abs(reference_point.y) || intersectionPoint.y <= reference_point.y)
			return -1;

		int tileX = (intersectionPoint.x - reference_point.x) / length;
		int tileY = (intersectionPoint.y - reference_point.y) / length;

		if (((tileX + tileY) % 2) == 0)
		{
			color[0] = 1;
			color[1] = 1;
			color[2] = 1;
		}
		else
		{
			color[0] = 0;
			color[1] = 0;
			color[2] = 0;
		}

		if (level == 0)
			return t;

		double intensity[3] = {0.4, 0.4, 0.4};

		double totalSpecularIntensity[3] = {0, 0, 0};
		for (int i = 0; i < lights.size(); i++)
		{
			// Diffuse
			Vector3D L = vectorSub(intersectionPoint, lights[i].light_pos);
			L = normalize(L);
			double diffuseIntensity = -0.2 * vectorDot(L, normal);

			// Specular
			Vector3D H = scale(vectorAdd(L, ray.dir), 0.5);
			H = normalize(H);
			double specDot = vectorDot(normal, H);
			int sign = 1;
			if (specDot < 0)
				sign = -1;
			specDot *= sign;
			double specularIntensity = -0.1 * sign * power(specDot, shine);

			for (int j = 0; j < 3; j++)
			{
				intensity[j] += max(0.0, diffuseIntensity) * lights[i].color[j];
				totalSpecularIntensity[j] += max(0.0, specularIntensity) * lights[i].color[j];
			}
		}

		for (int i = 0; i < 3; i++)
		{
			color[i] = min(1.0, color[i] * intensity[i] + totalSpecularIntensity[i]);
			color[i] = max(color[i], 0.0);
		}

		Vector3D reflectionDir = vectorSub(ray.dir, scale(normal, 2 * vectorDot(ray.dir, normal)));
		reflectionDir = normalize(reflectionDir);
		Ray reflectRay(vectorAdd(intersectionPoint, scale(reflectionDir, 0.2)), reflectionDir);
		//print(reflectionDir, "\n\nreflect dir");
		//print(intersectionPoint, "\n\nintersect point");
		//print(reference_point, "\n\nreference point");

		float minT = -1;
		int nearest = -1;
		for (int k = 0; k < objects.size(); k++)
		{
			if (this == objects[k])
				continue;
			double *dummyColor = new double[3];
			float newT = objects[k]->intersect(reflectRay, dummyColor, 0);
			if ((nearest == -1 || newT < minT) && newT > 0)
			{
				nearest = k;
				minT = newT;
			}
		}
		if (minT > 0)
		{
			double *nextColor = new double[3];
			minT = objects[nearest]->intersect(reflectRay, nextColor, level + 1);
			for (int i = 0; i < 3; i++)
			{
				//cout << nextColor[i] << " nextcolor \n";
				color[i] += nextColor[i] * 0.2;
				color[i] = min(color[i], 1.0);
				color[i] = max(color[i], 0.0);
			}
		}
		return t;
	}
};





void loadData(){
    ifstream sceneInput("scene.txt"); 
	string input;
	getline(sceneInput, input);

	int resolution;
	int objectCount;

	stringstream lineRecursion(input);
	lineRecursion >> recursionLevel;
	cout << "REC " << recursionLevel;

	getline(sceneInput, input);
	stringstream lineRes(input);
	lineRes >> resolution;

	imageWidth = resolution;
	imageHeight = resolution;


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
			lineRef >> object->length;
			lineRef >> object->width;
			lineRef >> object->height;

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

			//cout << "A" << object->A << endl;
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
	cout << "Light Count " << lightCount << endl;

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
	floor = new Floor(400, 20);
	objects.push_back(floor);

}
