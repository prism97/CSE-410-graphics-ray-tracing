//
// Created by ubuntu on ১৬/৩/২১.
//
//#include <windows.h>



#include "1605094_Header.h"
#include "bitmap_image.hpp"

using namespace std;



#define MOVE_CONST 2
#define ROTATE_CONST ((0.5/180) * PI) // 0.5 degree

int windowWidth = 500, windowHeight = 500;
double viewAngle = 80;
int recursion_level, pixels;
vector<Object*> objects;
vector<Light> lights;


int drawgrid;
int drawaxes;
//
//struct point
//{
//    double x,y,z;
//};





// camera variables
Vector3D eye, upVec, rightVec, lookVec;
//
//// utility functions
//void printVector(struct point A) {
//    printf("%lf, %lf, %lf\n", A.x, A.y, A.z);
//}
//
//struct point addVector(struct point A, struct point B) {
//    return {A.x + B.x, A.y + B.y, A.z + B.z};
//}
//
//struct point subtractVector(struct point A, struct point B) {
//    return {A.x - B.x, A.y - B.y, A.z - B.z};
//}
//
//void normalizeVector(struct point * A) {
//    double x = A -> x;
//    double y = A -> y;
//    double z = A -> z;
//    double magnitude = sqrt(x*x + y*y + z*z);
//    A -> x = x / magnitude;
//    A -> y = y / magnitude;
//    A -> z = z / magnitude;
//}
//
//struct point crossProduct(struct point A, struct point B) {
//    return {A.y * B.z - A.z * B.y,
//            A.z * B.x - A.x * B.z,
//            A.x * B.y - A.y * B.x};
//}
//
//struct point scalarMultiplication(double s, struct point A) {
//    return {s * A.x, s * A.y, s * A.z};
//}
//
//struct point rotateVector(struct point r, struct point v, double angle) {
//    // A rotates about B by an angle of theta
//    struct point res = addVector(scalarMultiplication(cos(angle), v),
//                                 scalarMultiplication(sin(angle), crossProduct(r, v)));
//    normalizeVector(&res);
//    return res;
//}



void drawAxes()
{
    if(drawaxes==1)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);{
            glVertex3f( 500,0,0);
            glVertex3f(-500,0,0);

            glVertex3f(0,500,0);
            glVertex3f(0, -500,0);

            glVertex3f(0,0,500);
            glVertex3f(0, 0,-500);
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


void drawSS()
{
    for (Object* o : objects) {
        o->draw();
    }
}

void display(){

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-upVec camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    gluLookAt(eye.x, eye.y, eye.z, eye.x + lookVec.x, eye.y + lookVec.y, eye.z + lookVec.z, upVec.x, upVec.y, upVec.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawGrid();

    drawSS();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}

void capture() {
    int imageWidth = pixels, imageHeight = pixels;
    // initialize bitmap image
    bitmap_image image(pixels,pixels);

    // set background color
    for (int i = 0; i < imageWidth; i++) {
        for (int j = 0; j < imageHeight; j++) {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }

    double planeDistance = (windowHeight / 2.0) / tan(viewAngle / 2.0);
    Vector3D topLeft = eye + (lookVec * planeDistance) - (rightVec * (windowWidth / 2.0)) + (upVec * (windowHeight / 2.0));
    double du = windowWidth / imageWidth;
    double dv = windowHeight / imageHeight;
    // choose middle of the grid cell
    topLeft = topLeft + (rightVec * (0.5 * du)) - (upVec * (0.5 * dv));

    int nearest;
    double t, tMin;
    Vector3D curPixel{};
    for (int i = 0; i < imageWidth; i++) {
        for (int j = 0; j < imageHeight; j++) {
            // calculate curPixel
            curPixel = topLeft + (rightVec * (i * du)) - (upVec * (j * dv));
            // cast ray
            Ray* ray = new Ray(eye, (curPixel - eye));
            double* color = new double[3];
            for (Object* o : objects) {
                t = o->intersect(ray, color, 0);

            }
            image.set_pixel(i, j, color[0], color[1], color[2]);
        }
    }
    image.save_image("output.bmp");;

}



void keyboardListener(unsigned char key, int x,int y){
    switch(key){
        case '0':
            capture();
            break;
        case '1':
            rightVec = rightVec.rotate(upVec, ROTATE_CONST);
            lookVec = lookVec.rotate(upVec, ROTATE_CONST);
            break;
        case '2':
            rightVec = rightVec.rotate(upVec, -ROTATE_CONST);
            lookVec = lookVec.rotate(upVec, -ROTATE_CONST);
            break;
        case '3':
            upVec = upVec.rotate(rightVec, ROTATE_CONST);
            lookVec = lookVec.rotate(rightVec, ROTATE_CONST);
            break;
        case '4':
            upVec = upVec.rotate(rightVec, -ROTATE_CONST);
            lookVec = lookVec.rotate(rightVec, -ROTATE_CONST);
            break;
        case '5':
            upVec = upVec.rotate(lookVec, -ROTATE_CONST);
            rightVec = rightVec.rotate(lookVec, -ROTATE_CONST);
            break;
        case '6':
            upVec = upVec.rotate(lookVec, ROTATE_CONST);
            rightVec = rightVec.rotate(lookVec, ROTATE_CONST);
            break;
        default:
            break;
    }
}

void specialKeyListener(int key, int x,int y){
    switch(key){
        case GLUT_KEY_UP:		// upVec arrow key
            eye = eye + (lookVec * MOVE_CONST);
            break;
        case GLUT_KEY_DOWN:		//down arrow key
            eye = eye - (lookVec * MOVE_CONST);
            break;
        case GLUT_KEY_RIGHT:
            eye = eye + (rightVec * MOVE_CONST);
            break;
        case GLUT_KEY_LEFT:
            eye = eye - (rightVec * MOVE_CONST);
            break;
        case GLUT_KEY_PAGE_UP:
            eye = eye + (upVec * MOVE_CONST);
            break;
        case GLUT_KEY_PAGE_DOWN:
            eye = eye - (upVec * MOVE_CONST);
            break;
        default:
            break;
    }
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
    switch(button){
        case GLUT_LEFT_BUTTON:
            if(state == GLUT_DOWN) {

            }
            break;
        case GLUT_RIGHT_BUTTON:
            if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
                drawaxes=1-drawaxes;
            }
            break;
        default:
            break;
    }
}




void animate(){
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}


void init() {
    //codes for initialization
    drawgrid=0;
    drawaxes=1;

    // initialize camera variables
    eye.x = 200;
    eye.y = 200;
    eye.z = 0;

    upVec.x = 0;
    upVec.y = 0;
    upVec.z = 1;
    upVec.normalize();

    rightVec.x = -1;
    rightVec.y = 1;
    rightVec.z = 0;
    rightVec.normalize();

    lookVec.x = -1;
    lookVec.y = -1;
    lookVec.z = 0;
    lookVec.normalize();

    //clear the screen
    glClearColor(0,0,0,0);

    /************************
    / set-upVec projection here
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
}

void loadData() {
    ifstream sceneFile;
    sceneFile.open("scene.txt");

    if (!sceneFile) {
        cerr << "Unable to open file scene.txt" << endl;
        exit(1);
    }

    int object_count;
    sceneFile >> recursion_level >> pixels >> object_count;

    string object_name;
    double x, y, z, w;
    Object *temp = nullptr;
    for (int i = 0; i < object_count; i++) {
        sceneFile >> object_name;
        if (object_name == "sphere") {
            sceneFile >> x >> y >> z;
            Vector3D center(x, y, z);
            double radius;
            sceneFile >> radius;
            temp = new Sphere(center, radius);
        } else if (object_name == "triangle") {
            Vector3D points[3];
            for (int j = 0; j < 3; j++) {
                sceneFile >> x >> y >> z;
                Vector3D point(x, y, z);
                points[j] = point;
            }
            temp = new Triangle(points);
        } else if (object_name == "general") {
            double A, B, C, D, E, F, G, H, I, J;
            sceneFile >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
            sceneFile >> x >> y >> z;
            Vector3D cube_ref(x, y, z);
            double length, width, height;
            sceneFile >> length >> width >> height;
            temp = new GeneralQuadricSurface();
        }
        sceneFile >> x >> y >> z;
        temp->setColor(x, y, z);
        sceneFile >> x >> y >> z >> w;
        temp->setCoefficients(x, y, z, w);
        int shine;
        sceneFile >> shine;
        temp->setShine(shine);

        objects.push_back(temp);
    }

    temp = new Floor(1000, 20);
    temp->setCoefficients(0.1, 0.1, 0.1, 0.1);
    temp->setShine(1);
    objects.push_back(temp);

    int light_count;
    sceneFile >> light_count;
    Light *l;
    for (int i = 0; i < light_count; i++) {
        sceneFile >> x >> y >> z;
        Vector3D position(x, y, z);
        double color[3];
        sceneFile >> color[0] >> color[1] >> color[2];
        l = new Light(position, color);
        lights.push_back(*l);
    }
}

int main(int argc, char **argv){
    loadData();

    glutInit(&argc,argv);
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occurring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}
