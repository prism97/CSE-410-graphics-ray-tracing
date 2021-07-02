#include <iostream>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <GL/glut.h>

using namespace std;

#define PI (2*acos(0.0))



class Vector3D {
public:
    double x, y, z;

    Vector3D()= default;
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
    Vector3D(const Vector3D &v) = default;

    Vector3D operator+(const Vector3D &v) const;
    Vector3D operator-(const Vector3D &v) const;
    Vector3D operator*(double s) const;
    Vector3D& operator=(const Vector3D &v);
    Vector3D rotate(const Vector3D& r, double angle);
    void normalize();
    void print();

    static double dotProduct(const Vector3D &a, const Vector3D &b);
    static Vector3D crossProduct(const Vector3D &a, const Vector3D &b);
};

Vector3D& Vector3D::operator=(const Vector3D &v) = default;

void Vector3D::normalize() {
    double len = sqrt(x*x + y*y + z*z);
    x = x / len;
    y = y / len;
    z = z / len;
}

void Vector3D::print() {
    cout << x << " " << y << " " << z << endl;
}

double Vector3D::dotProduct(const Vector3D &a, const Vector3D &b) {
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

Vector3D Vector3D::crossProduct(const Vector3D &a, const Vector3D &b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

Vector3D Vector3D::operator+(const Vector3D &v) const {
    return {x + v.x, y + v.y, z + v.z};
}

Vector3D Vector3D::operator-(const Vector3D &v) const {
    return {x - v.x, y - v.y, z - v.z};
}

Vector3D Vector3D::operator*(double s) const {
    return {s * x, s * y, s * z};
}

Vector3D Vector3D::rotate(const Vector3D &r, double angle) {
    Vector3D res = ((*this) * cos(angle)) + (crossProduct(r, *this) * sin(angle));
    res.normalize();
    return res;
}


class Ray {
    Vector3D start;
    Vector3D dir;

public:
    Ray(const Vector3D &start, const Vector3D &dir) {
        this->start = start;
        this->dir = dir;
        this->dir.normalize();
    }
};



class Object {
protected:
    Vector3D reference_point{};
    double height{}, width{}, length{};
    double color[3]{};
    double coefficients[4]{}; // reflection coefficients
    int shine{}; // exponent term of specular component

public:
    Object()= default;
    virtual void draw(){
    }
    virtual double intersect(Ray *r, double *color, int level){
        return -1.0;
    }
    void setColor(double c1, double c2, double c3){
        color[0] = c1;
        color[1] = c2;
        color[2] = c3;
    }
    void setShine(int s){
        shine = s;
    }
    void setCoefficients(double c1, double c2, double c3, double c4){
        coefficients[0] = c1;
        coefficients[1] = c2;
        coefficients[2] = c3;
        coefficients[3] = c4;
    }
    void print() {
        reference_point.print();
        cout << "color : " << color[0] << " " << color[1] << " " << color[2] << endl;
    }
};


class Light {
    Vector3D light_pos{};
    double color[3]{};

public:
    Light(const Vector3D &pos, const double c[3]) {
        light_pos = pos;
        for (int i = 0; i < 3; i++) {
            color[i] = c[i];
        }
    }
};

extern vector<Object*> objects;
extern vector<Light> lights;


class Sphere : public Object {
public:
    Sphere(Vector3D center, double radius) {
        reference_point = center;
        length = radius;
    }

    void draw() override {
        glTranslatef(reference_point.x, reference_point.y, reference_point.z);
        int stacks = 20, slices = 24;
        Vector3D points[100][100];
        int i,j;
        double h,r;
        double radius = length;
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(PI/2));
            r=radius*cos(((double)i/(double)stacks)*(PI/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*PI);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*PI);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            glColor3f(color[0],color[1],color[2]);
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
        glTranslatef(-reference_point.x, -reference_point.y, -reference_point.z);
    }

    double intersect(Ray *r, double *color, int level) override {
        return 0;
    }
};

class Triangle : public Object {
    Vector3D points[3]{};
public:
    explicit Triangle(Vector3D p[3]) {
        for (int i = 0; i < 3; i++) {
            points[i] = p[i];
        }
    }

    void draw() override {
        glBegin(GL_TRIANGLES);{
            glColor3f(color[0], color[1], color[2]);
            glVertex3f(points[0].x, points[0].y, points[0].z);
            glVertex3f(points[1].x, points[1].y, points[1].z);
            glVertex3f(points[2].x, points[2].y, points[2].z);
        }glEnd();
    }

    double intersect(Ray *r, double *color, int level) override {
        return 0;
    }
};

class GeneralQuadricSurface : public Object {
public:
    GeneralQuadricSurface() {}

    double intersect(Ray *r, double *color, int level) override {
        return 0;
    }
};


class Floor : public Object {
public:
    Floor(double floorWidth, double tileWidth) {
        reference_point = {- floorWidth / 2, - floorWidth / 2, 0};
        length = tileWidth;
    }

    void draw() override {
        glBegin(GL_QUADS);
        {
            int limit = - (int)(reference_point.x / length);
            for (int i = -limit; i < limit; i++) {
                for (int j = -limit; j < limit; j++) {
                    if ((i + j) % 2 == 0) glColor3f(1, 1, 1);
                    else glColor3f(0, 0, 0);
                    glVertex3f(i * length, j * length, 0);
                    glVertex3f(i * length + length, j * length, 0);
                    glVertex3f(i * length + length, j * length + length, 0);
                    glVertex3f(i * length, j * length + length, 0);
                }
            }
        }
        glEnd();
    }

    double intersect(Ray *r, double *color, int level) override {
        return 0;
    }
};





