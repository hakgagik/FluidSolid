#ifndef _IOSTREAM_
#include <iostream>
#define _IOSTREAM_
#endif

#ifndef _GLUT_
#include "Dependencies\glew\glew.h"
#include "Dependencies\freeglut\freeglut.h"
#define _GLUT_
#endif

#include "MACGrid.h"

#define PI 3.14159265359f

MACGrid grid;

bool simulate;

int windowWidth;
int windowHeight;
float viewR;
float viewTheta;
float viewPhi;
float dAngle = 0.1f;
float dScale = 1.1f;


void Simulate() {
	grid.timeStep();
}

void Display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	float x = viewR * sin(viewTheta) * cos(viewPhi);
	float y = viewR * sin(viewTheta) * sin(viewPhi);
	float z = viewR * cos(viewTheta);

	gluLookAt(x, y, z, (float)grid.length / 2.0, (float)grid.width / 2.0, (float)grid.height / 2.0, 0, 0, 1);

	grid.display();

	glutSwapBuffers();
}

void Reshape(int width, int height) {
	windowWidth = width;
	windowHeight = height;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(5, (float)windowWidth / (float)windowHeight, 0.01, viewR * 2);
	glMatrixMode(GL_MODELVIEW);
}

void Special(int key, int, int) {
	switch (key) {
	case GLUT_KEY_UP:
		viewTheta -= dAngle;
		if (viewTheta < 2 * dAngle) viewTheta = 2 * dAngle;
		break;
	case GLUT_KEY_DOWN:
		viewTheta += dAngle;
		if (viewTheta > PI - 2 * dAngle) viewTheta = PI - 2 * dAngle;
		break;
	case GLUT_KEY_RIGHT:
		viewPhi += dAngle;
		break;
	case GLUT_KEY_LEFT:
		viewPhi -= dAngle;
	}
}

void Keyboard(unsigned char key, int, int) {
	switch (key) {
	case '=':
		viewR /= dScale;
		break;
	case '-':
		viewR *= dScale;
		break;
	case ' ':
		simulate = !simulate;
		break;
	case 's':
		Simulate();
		glutPostRedisplay();
	}

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(5, (float)windowWidth / (float)windowHeight, 0.01, viewR * 2);
	glMatrixMode(GL_MODELVIEW);
}

void Timer(int t) {
	if (simulate) {
		Simulate();
	}
	glutPostRedisplay();
	glutTimerFunc(t, Timer, t);
}

void Menu(int t)
{
	
}

void init() {
	MACGrid::initGL("vert.glsl", "frag.glsl");
	grid = MACGrid::buildMacGrid("test.txt");
	//grid.printFlags();
	viewR = sqrt(grid.length * grid.length + grid.width + grid.width + grid.height * grid.height) * 20;
	viewTheta = (float)acos(1.0 / sqrt(3));
	viewPhi = -45.0 / 180.0 * PI;

	float lightPos[] = { (float)grid.length / 2.0f, (float)grid.width / 2.0f, (float)grid.height / 2.0f, 1.0f };
	float lightAmbient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
	float lightDiffuse[] = { 0.9f, 0.9f, 0.9f, 1.0f };
	float lightSpecular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
	glEnable(GL_LIGHT0);
	glEnable(GL_ALPHA);
	glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
}

int main(int argc, char** argv)
{	
	glutInit(&argc, argv);
	glutInitWindowSize(800, 800);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow("FluidSolid");
	if (glewInit() != GLEW_OK) throw std::runtime_error("glewInit failed");
	glutDisplayFunc(Display);
	glutReshapeFunc(Reshape);
	glutSpecialFunc(Special);
	glutKeyboardFunc(Keyboard);
	glutCreateMenu(Menu);
	glutAddMenuEntry("blah", 1);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	init();

	Timer(16);
	glutMainLoop();
}