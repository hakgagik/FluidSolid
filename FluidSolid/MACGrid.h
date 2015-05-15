#ifndef _IOSTREAM_
#include <iostream>
#define _IOSTREAM_
#endif

#ifndef _GLUT_
#include "Dependencies\glew\glew.h"
#include "Dependencies\freeglut\freeglut.h"
#define _GLUT_
#endif

class MACGrid {
public:
	static const float g;
	static const int particlesPerCell;

	static MACGrid buildMacGrid(std::string inFile);
	static MACGrid buildVortex(int, int, int, float, float);
	static void initGL(const char* vertPath, const char* fragPath);

	static bool doProject;
	static bool doAdvect;
	static bool showGrid;
	static bool showVel;
	static bool showParticles;
	static bool showPressure;

	int length;
	int width;
	int height;
	int elemsPerPlane;
	int numElements;
	int numParticles;

	float dt;
	float t;

	MACGrid();
	MACGrid(int, int, int);
	//~MACGrid();

	void timeStep();
	void display();
	void printFlags();
	void printU();
	void printV();
	void printW();

	inline float T2nu(float);
	inline float getU(float, float, float);
	inline float getV(float, float, float);
	inline float getW(float, float, float);
	inline float getUB(float, float, float);
	inline float getVB(float, float, float);
	inline float getWB(float, float, float);
	inline float interpolate(float***, float, float, float);
	
private:
	static const float particleRadius;
	static const int neighbors_I[8];
	static const int neighbors_J[8];
	static const int neighbors_K[8];
	static const int pressureProjectionSteps;
	static float relax;
	static int PARTICLE_DISPLAY_LIST;
	static GLuint program;

	float*** pressure;
	float*** T;
	float*** u;
	float*** v;
	float*** w;
	float*** uB;
	float*** vB;
	float*** wB;
	float*** nu;
	float*** delDotU;
	char*** flags;

	float* particlesX;
	float* particlesY;
	float* particlesZ;

	float maxSpeed;

	void getdt();
	void addforce();
	void advect();
	void diffuse();
	void project();
	void advectEverythingElse();
};