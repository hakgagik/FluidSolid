#include "MACGrid.h"
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include "Program.h"

#ifndef _GLUT_
#include "Dependencies\glew\glew.h"
#include "Dependencies\freeglut\freeglut.h"
#define _GLUT_
#endif

const int MACGrid::neighbors_I[8] = { 1, -1, 0, 0, 0, 0 };
const int MACGrid::neighbors_J[8] = { 0, 0, 1, -1, 0, 0 };
const int MACGrid::neighbors_K[8] = { 0, 0, 0, 0, 1, -1 };
const int MACGrid::pressureProjectionSteps = 1000;
const int MACGrid::particlesPerCell = 2;
const float MACGrid::g = -9.81f;
const float MACGrid::particleRadius = 0.3f;
float MACGrid::relax = 1.6f;
bool MACGrid::doProject = true;
bool MACGrid::doAdvect = true;
bool MACGrid::showGrid = false;
bool MACGrid::showVel = false;
bool MACGrid::showParticles = true;
bool MACGrid::showPressure = false;
int MACGrid::PARTICLE_DISPLAY_LIST = -1;
GLuint MACGrid::program = -1;

float deltaX, deltaY, deltaZ;
float dx, dy, dz;

MACGrid::MACGrid() {}

MACGrid::MACGrid(int length, int width, int height) {
	this->length = length;
	this->width = width;
	this->height = height;
	elemsPerPlane = length * width;
	numElements = height * elemsPerPlane;

	//Inintialize all the arrays
	pressure = new float**[length];
	flags = new char**[length];
	T = new float**[length];
	nu = new float**[length];
	delDotU = new float**[length];
	for (int i = 0; i < length; i++) {
		pressure[i] = new float*[width];
		flags[i] = new char*[width];
		T[i] = new float*[width];
		nu[i] = new float*[width];
		delDotU[i] = new float*[height];
		for (int j = 0; j < width; j++) {
			pressure[i][j] = new float[height];
			flags[i][j] = new char[height];
			T[i][j] = new float[height];
			nu[i][j] = new float[height];
			delDotU[i][j] = new float[height];
		}
	}

	u = new float**[length + 1];
	uB = new float**[length + 1];
	for (int i = 0; i <= length; i++) {
		u[i] = new float*[width];
		uB[i] = new float*[width];
		for (int j = 0; j < width; j++) {
			u[i][j] = new float[height];
			uB[i][j] = new float[height];
		}
	}	
	v = new float**[length];
	vB = new float**[length];
	for (int i = 0; i < length; i++) {
		v[i] = new float*[width + 1];
		vB[i] = new float*[width + 1];
		for (int j = 0; j <= width; j++) {
			v[i][j] = new float[height];
			vB[i][j] = new float[height];
		}
	}
	w = new float**[length];
	wB = new float**[length];
	for (int i = 0; i < length; i++) {
		w[i] = new float*[width];
		wB[i] = new float*[width];
		for (int j = 0; j < width; j++) {
			w[i][j] = new float[height + 1];
			wB[i][j] = new float[height + 1];
		}
	}
	t = 0;

	//Set everything to 0
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				pressure[i][j][k] = 0;
				T[i][j][k] = 293;
				nu[i][j][k] = T2nu(293);
				delDotU[i][j][k] = 0;
			}
		}
	}
	for (int i = 0; i <= length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k <= height; k++) {
				if (i == length && j == width || j == width && k == height || k == height && i == length) {
					continue;
				}
				if (i == length) {
					u[i][j][k] = 0;
					uB[i][j][k] = 0;
				}
				else if (j == length) {
					v[i][j][k] = 0;
					vB[i][j][k] = 0;
				}
				else if (k == height) {
					w[i][j][k] = 0;
					wB[i][j][k] = 0;
				}
				else {
					u[i][j][k] = 0;
					uB[i][j][k] = 0;
					v[i][j][k] = 0;
					vB[i][j][k] = 0;
					w[i][j][k] = 0;
					wB[i][j][k] = 0;
				}
			}
		}
	}
}

//MACGrid::~MACGrid() {
//	delete[] T;
//	delete[] pressure;
//	delete[] u;
//	delete[] v;
//	delete[] w;
//	delete[] uB;
//	delete[] vB;
//	delete[] wB;
//	delete[] nu;
//	delete[] flags;
//
//	delete[] particlesX;
//	delete[] particlesY;
//	delete[] particlesZ;
//}

// Builds a MAC grid given an input file. The first line of the file should specify the number of cells in each direction.
// Each subsequent line should be a triple of ints that specify which cells start with fluid in them.
MACGrid MACGrid::buildMacGrid(std::string infile) {
	using namespace std;
	ifstream fin;
	fin.open(infile);
	if (!fin.good()) {
		return MACGrid(0, 0, 0);
	}
	int x, y, z;
	string line;
	string buffer;
	getline(fin, line);
	stringstream ss(line);
	ss >> buffer;
	int length = stoi(buffer);
	ss >> buffer;
	int width = stoi(buffer);
	ss >> buffer;
	int height = stoi(buffer);
	ss >> buffer;
	float pressure = stof(buffer);
	ss >> buffer;
	float T = stof(buffer);

	MACGrid grid(length, width, height);

	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				grid.flags[i][j][k] = 'e';
			}
		}
	}

	int numFluidCells = 0;
	while (!fin.eof()) {
		getline(fin, line);
		ss = stringstream(line);
		ss >> buffer;
		x = stoi(buffer);
		ss >> buffer;
		y = stoi(buffer);
		ss >> buffer;
		z = stoi(buffer);
		grid.flags[x][y][z] = 'f';
		grid.pressure[x][y][z] = pressure;
		grid.T[x][y][z] = T;
		grid.nu[x][y][z] = grid.T2nu(T);
		numFluidCells++;
	}

	fin.close();

	grid.numParticles = numFluidCells * particlesPerCell;

	int X, Y, Z;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				if (grid.flags[i][j][k] == 'e') {
					for (int n = 0; n < 8; n++) {
						X = i + neighbors_I[n];
						Y = j + neighbors_J[n];
						Z = k + neighbors_K[n];
						if (X >= 0 && Y >= 0 && Z >= 0 && X < length && Y < width && Z < height) {
							if (grid.flags[X][Y][Z] == 'f') {
								grid.flags[i][j][k] = 's';
								break;
							}
						}
					}
				}
			}
		}
	}

	grid.particlesX = new float[grid.numParticles];
	grid.particlesY = new float[grid.numParticles];
	grid.particlesZ = new float[grid.numParticles];

	int count = 0;

	for (int i = 0; i < grid.length; i++) {
		for (int j = 0; j < grid.width; j++) {
			for (int k = 0; k < grid.height; k++) {
				if (grid.flags[i][j][k] == 'f') {
					for (int p = 0; p < particlesPerCell; p++) {
						grid.particlesX[count] = float(i) + float(rand()) / static_cast<float>(RAND_MAX);
						grid.particlesY[count] = float(j) + float(rand()) / static_cast<float>(RAND_MAX);
						grid.particlesZ[count] = float(k) + float(rand()) / static_cast<float>(RAND_MAX);
						grid.particlesX[i] = (grid.particlesX[i] > 0 ? grid.particlesX[i] : 0) < length ? grid.particlesX[i] : length - 0.0001f;
						grid.particlesY[i] = (grid.particlesY[i] > 0 ? grid.particlesY[i] : 0) < width ? grid.particlesY[i] : width - 0.0001f;
						grid.particlesZ[i] = (grid.particlesZ[i] > 0 ? grid.particlesZ[i] : 0) < height ? grid.particlesZ[i] : height - 0.0001f;
						count++;
					}
				}
			}
		}
	}
	//grid.printFlags();
	return grid;
}

MACGrid MACGrid::buildVortex(int length, int width, int height, float r, float gamma) {
	MACGrid grid(length, width, height);
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				if (k < 3) grid.flags[i][j][k] = 'f';
				else grid.flags[i][j][k] = 'e';
			}
		}
	}
	grid.numParticles = length * width * height * particlesPerCell;

	int X, Y, Z;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				if (grid.flags[i][j][k] == 'e') {
					for (int n = 0; n < 8; n++) {
						X = i + neighbors_I[n];
						Y = j + neighbors_J[n];
						Z = k + neighbors_K[n];
						if (X >= 0 && Y >= 0 && Z >= 0 && X < length && Y < width && Z < height) {
							if (grid.flags[X][Y][Z] == 'f') {
								grid.flags[i][j][k] = 's';
								break;
							}
						}
					}
				}
			}
		}
	}
	for (int i = 0; i <= length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				float x = i - float(length) / 2.0f;
				float y = j + 0.5 - float(width) / 2.0f;
				if (sqrt(x*x + y*y) < r) grid.u[i][j][k] = -y * gamma;
			}
		}
	}
	
	for (int i = 0; i < length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k < height; k++) {
				float x = i + 0.5 - float(length) / 2.0f + r / 12.0;
				float y = j - float(width) / 2.0f + r / 12.0;
				if (sqrt(x*x + y*y) < r) grid.v[i][j][k] = x * gamma;
			}
		}
	}

	for (int i = 0; i <= length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				float x = i - float(length) / 2.0f - r / 12.0f;
				float y = j + 0.5 - float(width) / 2.0 - r / 12.0f;
				if (sqrt(x*x + y*y) < r) grid.u[i][j][k] += -y * gamma;
			}
		}
	}

	for (int i = 0; i < length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k < height; k++) {
				float x = i + 0.5 - float(length) / 2.0f - r / 12.0f;
				float y = j - float(width) / 2.0f - r / 12.0f;
				if (sqrt(x*x + y*y) < r) grid.v[i][j][k] += x * gamma;
			}
		}
	}


	grid.particlesX = new float[grid.numParticles];
	grid.particlesY = new float[grid.numParticles];
	grid.particlesZ = new float[grid.numParticles];

	int count = 0;

	for (int i = 0; i < grid.length; i++) {
		for (int j = 0; j < grid.width; j++) {
			for (int k = 0; k < grid.height; k++) {
				if (grid.flags[i][j][k] == 'f') {
					for (int p = 0; p < particlesPerCell; p++) {
						grid.particlesX[count] = float(i) + float(rand()) / static_cast<float>(RAND_MAX);
						grid.particlesY[count] = float(j) + float(rand()) / static_cast<float>(RAND_MAX);
						grid.particlesZ[count] = float(k) + float(rand()) / static_cast<float>(RAND_MAX);
						grid.particlesX[i] = (grid.particlesX[i] > 0 ? grid.particlesX[i] : 0) < length ? grid.particlesX[i] : length - 0.0001f;
						grid.particlesY[i] = (grid.particlesY[i] > 0 ? grid.particlesY[i] : 0) < width ? grid.particlesY[i] : width - 0.0001f;
						grid.particlesZ[i] = (grid.particlesZ[i] > 0 ? grid.particlesZ[i] : 0) < height ? grid.particlesZ[i] : height - 0.0001f;
						count++;
					}
				}
			}
		}
	}
	//grid.printFlags();
	return grid;
}

// Initialize GL things for drawing the MACGrid (I.e. shaders and display lists)
void MACGrid::initGL(const char* vertPath, const char* fragPath) {
	using namespace std;

	program = Program::createProgram(vertPath, fragPath, "Lambertian");

	if (PARTICLE_DISPLAY_LIST < 0) {
		int displayListIndex = glGenLists(1);
		GLUquadric* quadric = gluNewQuadric();
		glNewList(displayListIndex, GL_COMPILE);
		gluSphere(quadric, particleRadius, 8, 4);
		glEndList();
		gluDeleteQuadric(quadric);
		cout << "MADE DISPLAY LIST " << displayListIndex << " : " << glIsList(displayListIndex) << endl;
		PARTICLE_DISPLAY_LIST = displayListIndex;
	}

	//TODO: add shader stuff;
}

void MACGrid::timeStep() {
	getdt();
	addforce();
	advect();
	diffuse();
	project();
	advectEverythingElse();
}

void MACGrid::getdt() {
	maxSpeed = 0;
	float speed;
	bool e0, e1;
	for (int i = 0; i <= length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				e0 = e1 = true;
				if (i > 0) {
					if (flags[i - 1][j][k] != 'e') e0 = false;
				}
				if (i < length) {
					if (flags[i][j][k] != 'e') e1 = false;
				}
				if (e0 && e1) continue;
				speed = abs(u[i][j][k]);
				if (speed > maxSpeed) maxSpeed = speed;
			}
		}
	}
	for (int i = 0; i < length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k < height; k++) {
				e0 = e1 = true;
				if (j > 0) {
					if (flags[i][j - 1][k] != 'e') e0 = false;
				}
				if (j < width) {
					if (flags[i][j][k] != 'e') e1 = false;
				}
				if (e0 && e1) continue;
				speed = abs(v[i][j][k]);
				if (speed > maxSpeed) maxSpeed = speed;
			}
		}
	}
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k <= height; k++) {
				e0 = e1 = true;
				if (k > 0) {
					if (flags[i][j][k - 1] != 'e') e0 = false;
				}
				if (k < height) {
					if (flags[i][j][k] != 'e') e1 = false;
				}
				if (e0 && e1) continue;
				speed = abs(w[i][j][k]);
				if (speed > maxSpeed) maxSpeed = speed;
			}
		}
	}

	if (maxSpeed <= 0) dt = 0.01;
	else dt = 0.5 / maxSpeed;
	dt = std::min(0.01f, dt);
	t += dt;
	std::cout << "time: " << t << " dt: " << dt << " ";
}

void MACGrid::addforce() {
	float dv = dt * g * height;
	bool e0, e1;
	for (int i = 0; i <= length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				// TODO: add other forces as necessary
				e0 = e1 = true;
				if (i > 0) {
					if (flags[i - 1][j][k] != 'e') e0 = false;
				}
				if (i < length) {
					if (flags[i][j][k] != 'e') e1 = false;
				}
				if (e0 && e1) {
					uB[i][j][k] = getU(i, j + 0.5, k + 0.5);
					continue;
				}
				uB[i][j][k] = u[i][j][k];
				if (i == 0) {
					uB[i][j][k] = uB[i][j][k] > 0 ? uB[i][j][k] : 0;
				}
				if (i == length) {
					uB[i][j][k] = uB[i][j][k] > 0 ? 0 : uB[i][j][k];
				}
			}
		}
	}
	for (int i = 0; i < length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k < height; k++) {
				// TODO: add other forces as necessary
				e0 = e1 = true;
				if (j > 0) {
					if (flags[i][j - 1][k] != 'e') e0 = false;
				}
				if (j < width) {
					if (flags[i][j][k] != 'e') e1 = false;
				}
				if (e0 && e1) {
					vB[i][j][k] = getV(i+0.5, j, k+0.5);
					continue;
				}
				vB[i][j][k] = v[i][j][k];
				if (j == 0) {
					vB[i][j][k] = vB[i][j][k] > 0 ? vB[i][j][k] : 0;
				}
				if (j == width) {
					vB[i][j][k] = vB[i][j][k] > 0 ? 0 : vB[i][j][k];
				}
			}
		}
	}
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k <= height; k++) {
				// TODO: add other forces as necessary
				e0 = e1 = true;
				if (k > 0) {
					if (flags[i][j][k - 1] != 'e') e0 = false;
				}
				if (k < height) {
					if (flags[i][j][k] != 'e') e1 = false;
				}
				if (e0 && e1) {
					wB[i][j][k] = getW(i + 0.5, j + 0.5, k);
					continue;
				}
				wB[i][j][k] = w[i][j][k] + dv;
				if (k == 0) {
					wB[i][j][k] = wB[i][j][k] > 0 ? wB[i][j][k] : 0;
				}
				if (k == height) {
					wB[i][j][k] = wB[i][j][k] > 0 ? 0 : wB[i][j][k];
				}
			}
		}
	}
}

void MACGrid::advect() {
	if (doAdvect) {
		float x, y, z;
		// u
		for (int i = 0; i <= length; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < height; k++) {
					x = float(i);
					y = j + 0.5f;
					z = k + 0.5f;
					u[i][j][k] = getUB(x - uB[i][j][k] * dt, y - getVB(x, y, z) * dt, z - getWB(x, y, z) * dt);
				}
			}
		}

		// v
		for (int i = 0; i < length; i++) {
			for (int j = 0; j <= width; j++) {
				for (int k = 0; k < height; k++) {
					x = i + 0.5f;
					y = float(j);
					z = k + 0.5f;
					v[i][j][k] = getVB(x - getUB(x, y, z) * dt, y - vB[i][j][k] * dt, z - getWB(x, y, z) * dt);
				}
			}
		}

		// w
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < height; k++) {
					x = i + 0.5f;
					y = j + 0.5f;
					z = float(k);
					w[i][j][k] = getWB(x - getUB(x, y, z) * dt, y - getVB(x, y, z) * dt, z - wB[i][j][k] * dt);
				}
			}
		}
	} else {
		for (int i = 0; i <= length; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < height; k++) {
					u[i][j][k] = uB[i][j][k];
				}
			}
		}

		// v
		for (int i = 0; i < length; i++) {
			for (int j = 0; j <= width; j++) {
				for (int k = 0; k < height; k++) {
					v[i][j][k] = vB[i][j][k];
				}
			}
		}

		// w
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < height; k++) {
					w[i][j][k] = wB[i][j][k];
				}
			}
		}
	}
}

void MACGrid::diffuse() {
	float n = 0;
	for (int i = 0; i <= length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				float count = 0;
				float accum = 0;
				if (i > 0) {
					count++;
					accum += u[i - 1][j][k];
				}
				if (i < length) {
					count++;
					accum += u[i + 1][j][k];
				}
				if (j > 0) {
					count++;
					accum += u[i][j - 1][k];
				}
				if (j < width - 1) {
					count++;
					accum += u[i][j + 1][k];
				}
				if (k > 0) {
					count++;
					accum += u[i][j][k - 1];
				}
				if (k < height -1 ) {
					count++;
					accum += u[i][j][k + 1];
				}
				accum -= count * u[i][j][k];
				accum *= n * dt;
				uB[i][j][k] = uB[i][j][k] + accum;
			}
		}
	}

	// v
	for (int i = 0; i < length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k < height; k++) {
				float count = 0;
				float accum = 0;
				if (i > 0) {
					count++;
					accum += v[i - 1][j][k];
				}
				if (i < length - 1) {
					count++;
					accum += v[i + 1][j][k];
				}
				if (j > 0) {
					count++;
					accum += v[i][j - 1][k];
				}
				if (j < width) {
					count++;
					accum += v[i][j + 1][k];
				}
				if (k > 0) {
					count++;
					accum += v[i][j][k - 1];
				}
				if (k < height - 1) {
					count++;
					accum += v[i][j][k + 1];
				}
				accum -= count * v[i][j][k];
				accum *= n * dt;
				vB[i][j][k] = v[i][j][k] + accum;
			}
		}
	}

	// w
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				float count = 0;
				float accum = 0;
				if (i > 0) {
					count++;
					accum += w[i - 1][j][k];
				}
				if (i < length - 1) {
					count++;
					accum += w[i + 1][j][k];
				}
				if (j > 0) {
					count++;
					accum += w[i][j - 1][k];
				}
				if (j < width - 1) {
					count++;
					accum += w[i][j + 1][k];
				}
				if (k > 0) {
					count++;
					accum += w[i][j][k - 1];
				}
				if (k < height) {
					count++;
					accum += w[i][j][k + 1];
				}
				accum -= count * w[i][j][k];
				accum *= n * dt;
				wB[i][j][k] = w[i][j][k] + accum;
			}
		}
	}
	for (int i = 0; i <= length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				u[i][j][k] = uB[i][j][k];
			}
		}
	}

	// v
	for (int i = 0; i < length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k < height; k++) {
				v[i][j][k] = vB[i][j][k];
			}
		}
	}

	// w
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				w[i][j][k] = wB[i][j][k];
			}
		}
	}

}

void MACGrid::project() {
	if (doProject) {
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < height; k++) {
					if (flags[i][j][k] == 'e') {
						delDotU[i][j][k] = 0;
					}
					else {
						delDotU[i][j][k] = u[i + 1][j][k] - u[i][j][k]
							+ v[i][j + 1][k] - v[i][j][k]
							+ w[i][j][k + 1] - w[i][j][k];
					}
				}
			}
		}

		float laplacian;
		float neighborCount;
		float maxDelta;
		float pOld;

		// Calculate pressure by SOR
		int step;
		for (step = 0; step < pressureProjectionSteps; step++) {
			maxDelta = 0;
			for (int i = 0; i < length; i++) {
				for (int j = 0; j < width; j++) {
					for (int k = 0; k < height; k++) {
						if (flags[i][j][k] != 'e') {
							neighborCount = 0;
							laplacian = 0;
							if (i > 0) {
								neighborCount++;
								laplacian += pressure[i - 1][j][k];
							}
							if (i + 1 < length) {
								neighborCount++;
								laplacian += pressure[i + 1][j][k];
							}
							if (j > 0) {
								neighborCount++;
								laplacian += pressure[i][j - 1][k];
							}
							if (j + 1 < width) {
								neighborCount++;
								laplacian += pressure[i][j + 1][k];
							}
							if (k > 0) {
								neighborCount++;
								laplacian += pressure[i][j][k - 1];
							}
							if (k + 1 < height) {
								neighborCount++;
								laplacian += pressure[i][j][k + 1];
							}
							if (neighborCount == 0) {
								pressure[i][j][k] = 0;
								continue;
							}
							laplacian -= delDotU[i][j][k];
							laplacian *= relax / neighborCount;
							pOld = pressure[i][j][k];
							pressure[i][j][k] = (1 - relax) * pressure[i][j][k] + laplacian;
							maxDelta = std::max(maxDelta, abs(pressure[i][j][k] - pOld));
						}
						else {
							pressure[i][j][k] = 0;
						}
					}
				}
			}
			if (maxDelta < 0.01) break;
		}
		std::cout << "steps: " << step << " maxDelta: " << maxDelta << std::endl;

		// Project incompressibility by subtracting u - grad(pressure)
		for (int i = 0; i <= length; i++) {
			for (int j = 0; j <= width; j++) {
				for (int k = 0; k <= height; k++) {
					if ((i == length) && (j == width) || (j == width && k == height) || (k == height && i == length)) {
						continue;
					}
					if (i == length) {
						// Only do u
						u[i][j][k] += pressure[i - 1][j][k];
						if (u[i][j][k] > 0) u[i][j][k] = 0;
					}
					else if (j == width) {
						// Only do v
						v[i][j][k] += pressure[i][j - 1][k];
						if (v[i][j][k] > 0) v[i][j][k] = 0;
					}
					else if (k == height){
						// Only do w
						w[i][j][k] += pressure[i][j][k - 1];
						if (w[i][j][k] > 0) w[i][j][k] = 0;
					}
					else {
						if (i == 0) {
							u[i][j][k] -= pressure[i][j][k];
							if (u[i][j][k] < 0) u[i][j][k] = 0;
						}
						else {
							u[i][j][k] -= pressure[i][j][k] - pressure[i - 1][j][k];
						}
						if (j == 0) {
							v[i][j][k] -= pressure[i][j][k];
							if (v[i][j][k] < 0) v[i][j][k] = 0;
						}
						else {
							v[i][j][k] -= pressure[i][j][k] - pressure[i][j - 1][k];
						}
						if (k == 0) {
							w[i][j][k] -= pressure[i][j][k];
							if (w[i][j][k] < 0) w[i][j][k] = 0;
						}
						else {
							w[i][j][k] -= pressure[i][j][k] - pressure[i][j][k - 1];
						}
					}
				}
			}
		}
	}
}

void MACGrid::advectEverythingElse() {
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				flags[i][j][k] = 'e';
			}
		}
	}
	for (int i = 0; i < numParticles; i++) {
		dx = getU(particlesX[i], particlesY[i], particlesZ[i]) * dt;
		dy = getV(particlesX[i], particlesY[i], particlesZ[i]) * dt;
		dz = getW(particlesX[i], particlesY[i], particlesZ[i]) * dt;
		particlesX[i] += dx;
		particlesY[i] += dy;
		particlesZ[i] += dz;
		particlesX[i] = particlesX[i] > 0 ? particlesX[i] : 0;
		particlesX[i] = particlesX[i] < length ? particlesX[i] : length - 0.0001f;
		particlesY[i] = particlesY[i] > 0 ? particlesY[i] : 0;
		particlesY[i] = particlesY[i] < width ? particlesY[i] : width - 0.0001f;
		particlesZ[i] = particlesZ[i] > 0 ? particlesZ[i] : 0;
		particlesZ[i] = particlesZ[i] < height ? particlesZ[i] : height - 0.0001f;
		flags[int(particlesX[i])][int(particlesY[i])][int(particlesZ[i])] = 'f';
	}
	// Check for surface cells. If a cell is a surface cell, the unallocated velocities around it become the average velocity of the allocated velocites around it.
	int X, Y, Z;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				if (flags[i][j][k] == 'e') {
					for (int n = 0; n < 8; n++) {
						X = i + neighbors_I[n];
						Y = j + neighbors_J[n];
						Z = k + neighbors_K[n];
						if (X >= 0 && Y >= 0 && Z >= 0 && X < length && Y < width && Z < height) {
							if (flags[X][Y][Z] == 'f') {
								flags[i][j][k] = 's';
							}
						}
					}
				}
			}
		}
	}
}

void MACGrid::display() {

	// Draw simulation box
	glColor3f(0, 0, 1);
	glBegin(GL_LINE_LOOP);
	glVertex3f(0, 0, 0);
	glVertex3f(0, width, 0);
	glVertex3f(length, width, 0);
	glVertex3f(length, 0, 0);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(0, 0, height);
	glVertex3f(0, width, height);
	glVertex3f(length, width, height);
	glVertex3f(length, 0, height);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, height);
	glVertex3f(0, width, 0);
	glVertex3f(0, width, height);
	glVertex3f(length, width, 0);
	glVertex3f(length, width, height);
	glVertex3f(length, 0, 0);
	glVertex3f(length, 0, height);
	glEnd();

	if (showPressure) {
		float max = -1E10;
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < height; k++) {
					if (pressure[i][j][k] > max) {
						max = pressure[i][j][k];
					}
				}
			}
		}
		if (max == 0) max = 1;
		glUseProgram(program);
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < height; k++) {
					float d[] = { pressure[i][j][k] / max, 0, 1 - pressure[i][j][k] / max, 1 };
					if (flags[i][j][k] == 'e') continue;
					glMaterialfv(GL_FRONT, GL_DIFFUSE, d);
					glPushMatrix();
					glTranslated(i + 0.5, j + 0.5, k + 0.5);
					glCallList(PARTICLE_DISPLAY_LIST);
					glPopMatrix();
				}
			}
		}
		glUseProgram(0);
	}

	if (showParticles) {
		glUseProgram(program);
		float c[] = { 0, 0, 1, 1 };
		for (int i = 0; i < numParticles; i++) {
			c[1] =  sqrt(particlesZ[i] / float(height));
			glMaterialfv(GL_FRONT, GL_DIFFUSE, c);
			glPushMatrix();
			glTranslated(particlesX[i], particlesY[i], particlesZ[i]);
			glCallList(PARTICLE_DISPLAY_LIST);
			glPopMatrix();
		}
		glUseProgram(0);
	}

	if (showVel) {
		float x, y, z, u, v, w;
		for (auto i = 0; i < length; i++) {
			for (auto j = 0; j < width; j++) {
				for (auto k = 0; k < height; k++) {
					if (flags[i][j][k] == 'e') continue;
					x = i + 0.5f;
					y = j + 0.5f;
					z = k + 0.5f;
					u = getU(x, y, z) / 10.0f;
					v = getV(x, y, z) / 10.0f;
					w = getW(x, y, z) / 10.0f;
					glBegin(GL_LINES);
					glColor3f(0, 0, 0);
					glVertex3d(x - u * 0.5, y - v * 0.5, z - w * 0.5);
					glColor3f(1, 1, 1);
					glVertex3d(x + u * 0.5, y + v * 0.5, z + w * 0.5);
					glEnd();
				}
			}
		}
	}

	if (showGrid) {
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < height; k++) {
					if (flags[i][j][k] != 'e') {
						if (flags[i][j][k] == 's') {
							glColor3f(1, 0, 0);
						}
						else if (flags[i][j][k] == 'f') {
							glColor3f(0, 1, 0);
						}
						else {
							glColor3f(0, 0, 0);
						}
						glBegin(GL_LINES);
						glVertex3f(i, j, k);
						if (i == length) {
							glVertex3f(i - 1, j, k);
						}
						else {
							glVertex3f(i + 1, j, k);
						}
						glVertex3f(i, j, k);
						if (j == width) {
							glVertex3f(i, j - 1, k);
						}
						else {
							glVertex3f(i, j + 1, k);
						}
						glVertex3f(i, j, k);
						if (k == height) {
							glVertex3f(i, j, k - 1);
						}
						else {
							glVertex3f(i, j, k + 1);
						}
						glEnd();
					}
				}
			}
		}
	}
}

void MACGrid::printFlags() {
	for (int k = 0; k < height; k++) {
		for (int j = 0; j < width; j++) {
			for (int i = 0; i < length; i++) {
				std::cout << flags[i][j][k];
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

}

void MACGrid::printU() {
	using namespace std;
	for (int k = 0; k < height; k++) {
		for (int j = 0; j < width; j++) {
			for (int i = 0; i <= length; i++) {
				cout << u[i][j][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

void MACGrid::printV() {
	using namespace std;
	for (int k = 0; k < height; k++) {
		for (int j = 0; j <= width; j++) {
			for (int i = 0; i < length; i++) {
				cout << v[i][j][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

void MACGrid::printW() {
	using namespace std;
	for (int k = 0; k <= height; k++) {
		for (int j = 0; j < width; j++) {
			for (int i = 0; i < length; i++) {
				cout << w[i][j][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

inline float MACGrid::T2nu(float T) {
	return 100.0f / T;
}

inline float MACGrid::getU(float x, float y, float z) {
	if (x < 0) x = 0;
	else if (x >= length) x = float(length) - 0.001f;
	if (y < 0.5) y = 0.5;
	else if (y >= width - 0.5) y = width - 0.5001f;
	if (z < 0.5) z = 0.5;
	else if (z >= height - 0.5) z = height - 0.5001f;
	int i = int(x);
	int j = int(y - 0.5f);
	int k = int(z - 0.5f);
	deltaX = x - int(x);
	deltaY = y - j - 0.5f;
	deltaZ = z - k - 0.5f;
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * u[i][j][k] + deltaX * u[i + 1][j][k])
		+ deltaY * ((1 - deltaX) * u[i][j + 1][k] + deltaX * u[i + 1][j + 1][k]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * u[i][j][k + 1] + deltaX * u[i + 1][j][k + 1])
		+ deltaY * ((1 - deltaX) * u[i][j + 1][k + 1] + deltaX * u[i + 1][j + 1][k + 1]));
}

inline float MACGrid::getV(float x, float y, float z) {
	if (x < 0.5) x = 0.5;
	else if (x >= length - 0.5) x = length - 0.5001f;
	if (y < 0) y = 0;
	else if (y >= width) y = float(width) - 0.001f;
	if (z < 0.5) z = 0.5;
	else if (z >= height - 0.5) z = height - 0.5001f;
	int i = int(x - 0.5f);
	int j = int(y);
	int k = int(z - 0.5f);
	deltaX = x - i - 0.5f;
	deltaY = y - j;
	deltaZ = z - k - 0.5f;
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * v[i][j][k] + deltaX * v[i + 1][j][k])
		+ deltaY * ((1 - deltaX) * v[i][j + 1][k] + deltaX * v[i + 1][j + 1 ][k]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * v[i][j][k + 1] + deltaX * v[i + 1][j][k])
		+ deltaY * ((1 - deltaX) * v[i][j + 1][k + 1] + deltaX * v[i + 1][j + 1][k + 1]));
}

inline float MACGrid::getW(float x, float y, float z) {
	if (x < 0.5) x = 0.5;
	else if (x >= length-0.5) x = length - 0.5001f;
	if (y < 0.5) y = 0.5;
	else if (y >= width - 0.5) y = width - 0.5001f;
	if (z < 0) z = 0;
	else if (z >= height) z = float(height) - 0.001f;
	int i = int(x - 0.5f);
	int j = int(y - 0.5f);
	int k = int(z);
	deltaX = x - i - 0.5f;
	deltaY = y - j - 0.5f;
	deltaZ = z - k;
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * w[i][j][k] + deltaX * w[i + 1][j][k])
		+ deltaY * ((1 - deltaX) * w[i][j + 1][k] + deltaX * w[i + 1][j + 1][k]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * w[i][j][k + 1] + deltaX * w[i + 1][j][k + 1])
		+ deltaY * ((1 - deltaX) * w[i][j + 1][k + 1] + deltaX * w[i + 1][j + 1][k + 1]));
}

inline float MACGrid::getUB(float x, float y, float z) {
	if (x < 0) x = 0;
	else if (x >= length) x = float(length) - 0.001f;
	if (y < 0.5) y = 0.5;
	else if (y >= width - 0.5) y = width - 0.5001f;
	if (z < 0.5) z = 0.5;
	else if (z >= height - 0.5) z = height - 0.5001f;
	int i = int(x);
	int j = int(y - 0.5f);
	int k = int(z - 0.5f);
	deltaX = x - int(x);
	deltaY = y - j - 0.5f;
	deltaZ = z - k - 0.5f;
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * uB[i][j][k] + deltaX * uB[i + 1][j][k])
		+ deltaY * ((1 - deltaX) * uB[i][j + 1][k] + deltaX * uB[i + 1][j + 1][k]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * uB[i][j][k + 1] + deltaX * uB[i + 1][j][k + 1])
		+ deltaY * ((1 - deltaX) * uB[i][j + 1][k + 1] + deltaX * uB[i + 1][j + 1][k + 1]));
}

inline float MACGrid::getVB(float x, float y, float z) {
	if (x < 0.5) x = 0.5;
	else if (x >= length - 0.5) x = length - 0.5001f;
	if (y < 0) y = 0;
	else if (y >= width) y = float(width) - 0.001f;
	if (z < 0.5) z = 0.5;
	else if (z >= height - 0.5) z = height - 0.5001f;
	int i = int(x - 0.5f);
	int j = int(y);
	int k = int(z - 0.5f);
	deltaX = x - i - 0.5f;
	deltaY = y - j;
	deltaZ = z - k - 0.5f;
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * vB[i][j][k] + deltaX * vB[i + 1][j][k])
		+ deltaY * ((1 - deltaX) * vB[i][j + 1][k] + deltaX * vB[i + 1][j + 1][k]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * vB[i][j][k + 1] + deltaX * vB[i + 1][j][k + 1])
		+ deltaY * ((1 - deltaX) * vB[i][j + 1][k + 1] + deltaX * vB[i + 1][j + 1][k + 1]));
}

inline float MACGrid::getWB(float x, float y, float z) {
	if (x < 0.5) x = 0.5;
	else if (x >= length - 0.5) x = length - 0.5001f;
	if (y < 0.5) y = 0.5;
	else if (y >= width - 0.5) y = width - 0.5001f;
	if (z < 0) z = 0;
	else if (z >= height) z = float(height) - 0.001f;
	int i = int(x - 0.5f);
	int j = int(y - 0.5f);
	int k = int(z);
	deltaX = x - i - 0.5f;
	deltaY = y - j - 0.5f;
	deltaZ = z - k;
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * wB[i][j][k] + deltaX * wB[i + 1][j][k])
		+ deltaY * ((1 - deltaX) * wB[i][j + 1][k] + deltaX * wB[i + 1][j + 1][k]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * wB[i][j][k + 1] + deltaX * wB[i + 1][j][k + 1])
		+ deltaY * ((1 - deltaX) * wB[i][j + 1][k + 1] + deltaX * wB[i + 1][j + 1][k + 1]));
}

inline float MACGrid::interpolate(float*** q, float x, float y, float z) {
	if (x < 0) x = 0;
	else if (x >= length - 1) x = length - 1.0001f;
	if (y < 0) y = 0;
	else if (y >= width - 1) y = width - 1.0001f;
	if (z < 0) z = 0;
	else if (z >= height - 1) z = height - 1.0001f;
	int i = int(x);
	int j = int(y);
	int k = int(z);
	deltaX = x - i;
	deltaY = y - j;
	deltaZ = z - k;
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * q[i][j][k] + deltaX * q[i + 1][j][k])
		+ deltaY * ((1 - deltaX) * q[i][j + 1][k] + deltaX * q[i + 1][j + 1][k]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * q[i][j][k + 1] + deltaX * q[i + 1][j][k + 1])
		+ deltaY * ((1 - deltaX) * q[i][j + 1][k + 1] + deltaX * q[i + 1][j + 1][k + 1]));
}