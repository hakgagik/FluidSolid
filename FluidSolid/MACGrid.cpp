#include "MACGrid.h"

#ifndef _PARSING_
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#define _PARSING_
#endif
#include <algorithm>

const int MACGrid::neighbors_I[8] = { 1, -1, 0, 0, 0, 0 };
const int MACGrid::neighbors_J[8] = { 0, 0, 1, -1, 0, 0 };
const int MACGrid::neighbors_K[8] = { 0, 0, 0, 0, 1, -1 };
const int MACGrid::pressureProjectionSteps = 50;
const int MACGrid::particlesPerCell = 2;
const float MACGrid::g = -9.81f;
const float MACGrid::particleRadius = 0.3f;
float MACGrid::relax = 1.6f;
bool MACGrid::debug = false;
bool MACGrid::showGrid = false;
int MACGrid::PARTICLE_DISPLAY_LIST = -1;
GLuint MACGrid::program = -1;

int ijk;
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
	pressure = new float[length * width * height];
	flags = new char[length * width * height];
	T = new float[length * width * height];
	u = new float[(length + 1) * width * height];
	v = new float[length * (width + 1) * height];
	w = new float[length * width * (height + 1)];
	uB = new float[(length + 1) * width * height];
	vB = new float[length * (width + 1) * height];
	wB = new float[length * width * (height + 1)];
	nu = new float[length * width * height];
	delDotU = new float[length * width * height];
	t = 0;

	//Set everything to 0
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				ijk = t1D0(i, j, k);
				pressure[ijk] = 0;
				T[ijk] = 293;
				nu[ijk] = T2nu(293);
				delDotU[ijk] = 0;
			}
		}
	}
	for (int i = 0; i <= length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k <= height; k++) {
				int ijk;
				if (i == length && j == width || j == width && k == height || k == height && i == length) {
					continue;
				}
				else if (i == length) {
					ijk = t1DU(i, j, k);
					u[ijk] = 0;
					uB[ijk] = 0;
				}
				else if (j == length) {
					ijk = t1DV(i, j, k);
					v[ijk] = 0;
					vB[ijk] = 0;
				}
				else if (k == height) {
					ijk = t1DW(i, j, k);
					w[ijk] = 0;
					wB[ijk] = 0;
				}
				else {
					ijk = t1DU(i, j, k);
					u[ijk] = 0;
					uB[ijk] = 0;
					ijk = t1DV(i, j, k);
					v[ijk] = 0;
					vB[ijk] = 0;
					ijk = t1DW(i, j, k);
					w[ijk] = 0;
					wB[ijk] = 0;
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

// Brings a MAC grid given an input file. The first line of the file should specify the number of cells in each direction.
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
				grid.flags[grid.t1D0(i, j, k)] = 'e';
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
		grid.flags[grid.t1D0(x, y, z)] = 'f';
		grid.pressure[grid.t1D0(x, y, z)] = pressure;
		grid.T[grid.t1D0(x, y, z)] = T;
		grid.nu[grid.t1D0(x, y, z)] = grid.T2nu(T);
		numFluidCells++;
	}

	fin.close();

	grid.numParticles = numFluidCells * particlesPerCell;

	int X, Y, Z;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				if (grid.flags[grid.t1D0(i, j, k)] == 'e') {
					for (int n = 0; n < 8; n++) {
						X = i + neighbors_I[n];
						Y = j + neighbors_J[n];
						Z = k + neighbors_K[n];
						if (X >= 0 && Y >= 0 && Z >= 0 && X < length && Y < width && Z < height) {
							if (grid.flags[grid.t1D0(X, Y, Z)] == 'f') {
								grid.flags[grid.t1D0(i, j, k)] = 's';
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
				int ijk = grid.t1D0(i, j, k);
				if (grid.flags[ijk] == 'f') {
					for (int p = 0; p < particlesPerCell; p++) {
						grid.particlesX[count] = (float)i + (float)rand() / static_cast<float>(RAND_MAX);
						grid.particlesY[count] = (float)j + (float)rand() / static_cast<float>(RAND_MAX);
						grid.particlesZ[count] = (float)k + (float)rand() / static_cast<float>(RAND_MAX);
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

std::string readShaderFile(const char *filePath) {
	using namespace std;
	string content;
	ifstream fileStream(filePath, ios::in);

	if (!fileStream.is_open()) {
		cerr << "Could not read file " << filePath << ". File does not exist." << endl;
		return "";
	}
	string line = "";
	while (!fileStream.eof()) {
		getline(fileStream, line);
		content.append(line + "\n");
	}
	fileStream.close();
	return content;
}

// Initialize GL things for drawing the MACGrid (I.e. shaders and display lists)
void MACGrid::initGL(const char* vertPath, const char* fragPath) {
	using namespace std;

	GLuint vertShader = glCreateShader(GL_VERTEX_SHADER);
	GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);

	// Read shader
	string vertShaderStr = readShaderFile(vertPath);
	string fragShaderStr = readShaderFile(fragPath);
	const char* vertShaderSrc = vertShaderStr.c_str();
	const char* fragShaderSrc = fragShaderStr.c_str();

	GLint result = GL_FALSE;
	int logLength;

	// Compile vertex shader
	cout << "Compiling vertex shader." << endl;
	glShaderSource(vertShader, 1, &vertShaderSrc, NULL);
	glCompileShader(vertShader);

	// Check vertex shader
	glGetShaderiv(vertShader, GL_COMPILE_STATUS, &result);
	glGetShaderiv(vertShader, GL_INFO_LOG_LENGTH, &logLength);
	vector<char> vertShaderError((logLength > 1) ? logLength : 1);
	glGetShaderInfoLog(vertShader, logLength, NULL, &vertShaderError[0]);
	cout << &vertShaderError[0] << std::endl;

	// Compile fragment shader
	cout << "Compiling fragment shader" << endl;
	glShaderSource(fragShader, 1, &fragShaderSrc, NULL);
	glCompileShader(fragShader);

	// Check fragment shader
	glGetShaderiv(fragShader, GL_COMPILE_STATUS, &result);
	glGetShaderiv(fragShader, GL_INFO_LOG_LENGTH, &logLength);
	vector<char> fragShaderError((logLength > 1) ? logLength : 1);
	glGetShaderInfoLog(fragShader, logLength, NULL, &fragShaderError[0]);
	std::cout << &fragShaderError[0] << std::endl;

	cout << "Linking program" << endl;
	program = glCreateProgram();
	glAttachShader(program, vertShader);
	glAttachShader(program, fragShader);
	glLinkProgram(program);

	glGetProgramiv(program, GL_LINK_STATUS, &result);
	glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);
	vector<char> programError((logLength > 1) ? logLength : 1);
	glGetProgramInfoLog(program, logLength, NULL, &programError[0]);
	cout << &programError[0] << std::endl;

	glDeleteShader(vertShader);
	glDeleteShader(fragShader);

	if (PARTICLE_DISPLAY_LIST < 0) {
		int displayListIndex = glGenLists(1);
		GLUquadric* quadric = gluNewQuadric();
		glNewList(displayListIndex, GL_COMPILE);
		gluSphere(quadric, particleRadius, 16, 8);
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
	//diffuse();
	project();
	advectEverythingElse();
}

void MACGrid::getdt() {
	float maxSpeed = 0;
	float speed;
	for (int i = 0; i <= length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				ijk = t1DU(i, j, k);
				speed = abs(u[ijk]);
				if (speed > maxSpeed) maxSpeed = speed;
			}
		}
	}
	for (int i = 0; i < length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k < height; k++) {
				ijk = t1DV(i, j, k);
				speed = abs(v[ijk]);
				if (speed > maxSpeed) maxSpeed = speed;
			}
		}
	}
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k <= height; k++) {
				ijk = t1DW(i, j, k);
				speed = abs(w[ijk]);
				if (speed > maxSpeed) maxSpeed = speed;
			}
		}
	}

	if (maxSpeed <= 0) dt = 0.01;
	else dt = 0.5 / maxSpeed;
	dt = std::min(0.2f, dt);
}

void MACGrid::addforce() {
	float dv = dt * MACGrid::g;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k <= height; k++) {
				// TODO: add other forces as necessary
				ijk = t1DW(i, j, k);
				w[ijk] += dv;
				if (k == 0) {
					w[ijk] = (w[ijk] > 0) ? w[ijk] : 0;
				}
			}
		}
	}
}

void MACGrid::advect() {
	float x, y, z;
	// u
	for (int i = 0; i <= length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				x = (float)i;
				y = j + 0.5f;
				z = k + 0.5f;
				uB[t1DU(i, j, k)] = getU(x - u[t1DU(i, j, k)] * dt, y - getV(x, y, z) * dt, z - getW(x, y, z) * dt);
			}
		}
	}

	// v
	for (int i = 0; i < length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k < height; k++) {
				x = i + 0.5f;
				y = (float)j;
				z = k + 0.5f;
				vB[t1DV(i, j, k)] = getV(x - getU(x, y, z) * dt, y - v[t1DV(i, j, k)] * dt, z - getW(x, y, z) * dt);
			}
		}
	}

	// w
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				x = i + 0.5f;
				y = j + 0.5f;
				z = (float)k;
				wB[t1DW(i, j, k)] = getW(x - getU(x, y, z) * dt, y - getV(x, y, z) * dt, z - w[t1DW(i, j, k)] * dt);
			}
		}
	}
}

void MACGrid::project() {
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				if (flags[t1D0(i, j, k)] == 'e') {
					delDotU[t1D0(i, j, k)] = 0;
				}
				else {
					delDotU[t1D0(i, j, k)] = uB[t1DU(i + 1, j, k)] - uB[t1DU(i, j, k)]
						+ vB[t1DV(i, j + 1, k)] - vB[t1DV(i, j, k)]
						+ wB[t1DW(i, j, k + 1)] - wB[t1DW(i, j, k)];
				}
			}
		}
	}

	float laplacian;
	float neighborCount;
	float maxDelta;
	float pOld;
	//std::cout << "Step" << std::endl;
	// Calculate pressure by SOR
	for (int step = 0; step < pressureProjectionSteps; step++) {
		maxDelta = 0;
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < height; k++) {
					ijk = t1D0(i, j, k);
					if (flags[ijk] != 'e') {
						neighborCount = 0;
						laplacian = 0;
						if (i - 1 >= 0) {
							neighborCount++;
							laplacian += pressure[ijk - 1];
						}
						if (i + 1 < length) {
							neighborCount++;
							laplacian += pressure[ijk + 1];
						}
						if (j - 1 >= 0) {
							neighborCount++;
							laplacian += pressure[ijk - length];
						}
						if (j + 1 < width) {
							neighborCount++;
							laplacian += pressure[ijk + length];
						}
						if (k - 1 >= 0) {
							neighborCount++;
							laplacian += pressure[ijk - elemsPerPlane];
						}
						if (k + 1 < height) {
							neighborCount++;
							laplacian += pressure[ijk + elemsPerPlane];
						}
						if (neighborCount == 0) {
							pressure[ijk] = 0;
							continue;
						}
						laplacian -= delDotU[ijk];
						laplacian *= relax / neighborCount;
						pOld = pressure[ijk];
						pressure[ijk] = (1 - relax) * pressure[ijk] + laplacian;
						maxDelta = std::max(maxDelta, abs(pressure[ijk] - pOld));
					}
					else {
						pressure[ijk] = 0;
					}
				}
			}
		}
		//std::cout << maxDelta << std::endl;
		if (maxDelta < 0.01) break;
	}

	// Project incompressibility by subtracting u - grad(pressure)
	for (int i = 0; i <= length; i++) {
		for (int j = 0; j <= width; j++) {
			for (int k = 0; k <= height; k++) {
				if ((i == length) && (j == width) || (j == width && k == height) || (k == height && i == length)) {
					continue;
				}
				if (i == length) {
					// Only do u
					u[t1DU(i, j, k)] -= pressure[t1D0(i - 1, j, k)];
					if (u[t1DU(i, j, k)] > 0) u[t1DU(i, j, k)] = 0;
				}
				else if (j == width) {
					// Only do v
					v[t1DV(i, j, k)] -= pressure[t1D0(i, j - 1, k)];
					if (v[t1DV(i, j, k)] > 0) v[t1DV(i, j, k)] = 0;
				}
				else if (k == height){
					// Only do w
					w[t1DW(i, j, k)] -= pressure[t1D0(i, j, k - 1)];
					if (w[t1DW(i, j, k)] > 0) w[t1DW(i, j, k)] = 0;
				}
				else {
					ijk = t1D0(i, j, k);
					if (i == 0) {
						u[t1DU(i, j, k)] -= pressure[ijk];
						if (u[t1DU(i, j, k)] < 0) u[t1DU(i, j, k)] = 0;
					}
					else {
						u[t1DU(i, j, k)] -= pressure[ijk] - pressure[ijk - 1];
					}
					if (j == 0) {
						v[t1DV(i, j, k)] -= pressure[ijk];
						if (v[t1DV(i, j, k)] < 0) v[t1DV(i, j, k)] = 0;
					}
					else {
						v[t1DV(i, j, k)] -= pressure[ijk] - pressure[ijk - length];
					}
					if (k == 0) {
						w[t1DW(i, j, k)] -= pressure[ijk];
						if (w[t1DW(i, j, k)] < 0) w[t1DW(i, j, k)] = 0;
					}
					else {
						w[t1DW(i, j, k)] -= pressure[ijk] - pressure[ijk - elemsPerPlane];
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
				ijk = t1D0(i, j, k);
				flags[ijk] = 'e';
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
		flags[t1D0((int)particlesX[i], (int)particlesY[i], (int)particlesZ[i])] = 'f';
	}
	int X, Y, Z;
	for (int i = 0; i < length; i++) {
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < height; k++) {
				if (flags[t1D0(i, j, k)] == 'e') {
					for (int n = 0; n < 8; n++) {
						X = i + neighbors_I[n];
						Y = j + neighbors_J[n];
						Z = k + neighbors_K[n];
						if (X >= 0 && Y >= 0 && Z >= 0 && X < length && Y < width && Z < height) {
							if (flags[t1D0(X, Y, Z)] == 'f') {
								flags[t1D0(i, j, k)] = 's';
								break;
							}
						}
					}
				}
			}
		}
	}
}

void MACGrid::display() {
	glUseProgram(program);
	float c[] = { 0, 0, 1, 1 };
	glMaterialfv(GL_FRONT, GL_DIFFUSE, c);
	for (int i = 0; i < numParticles; i++) {
		glPushMatrix();
		glTranslated(particlesX[i], particlesY[i], particlesZ[i]);
		glCallList(PARTICLE_DISPLAY_LIST);
		glPopMatrix();
	}
	glUseProgram(0);

	if (showGrid) {
		for (int i = 0; i < length; i++) {
			for (int j = 0; j < width; j++) {
				for (int k = 0; k < height; k++) {
					if (flags[t1D0(i, j, k)] != 'e') {
						if (flags[t1D0(i, j, k)] == 's') {
							glColor3f(1, 0, 0);
						}
						else if (flags[t1D0(i, j, k)] == 'f') {
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
				std::cout << flags[t1D0(i,j,k)];
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
				int ijk = t1DU(i, j, k);
				cout << u[ijk] << " ";
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
				int ijk = t1DV(i, j, k);
				cout << v[ijk] << " ";
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
				int ijk = t1DW(i, j, k);
				cout << w[ijk] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
}

inline int MACGrid::t1D0(int i, int j, int k) {
	return k * elemsPerPlane + j * length + i;
}

inline int MACGrid::t1DU(int i, int j, int k) {
	return k * (elemsPerPlane + width) + j * (length + 1) + i;
}

inline int MACGrid::t1DV(int i, int j, int k) {
	return k * (elemsPerPlane + length) + j * length + i;
}

inline int MACGrid::t1DW(int i, int j, int k) {
	return k * elemsPerPlane + j * length + i;
}

inline float MACGrid::T2nu(float T) {
	return 100.0f / T;
}

inline float MACGrid::getU(float x, float y, float z) {
	if (x < 0) x = 0;
	else if (x > length) x = (float)length;
	if (y < 0.5) y = 0.5;
	else if (y > width - 0.5) y = width - 0.5f;
	if (z < 0.5) z = 0.5;
	else if (z > height - 0.5) z = height - 0.5f;
	int i = (int)x;
	int j = (int)(y - 0.5f);
	int k = (int)(z - 0.5f);
	deltaX = x - (int)x;
	deltaY = y - j - 0.5f;
	deltaZ = z - k - 0.5f;
	ijk = t1DU(i, j, k);
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * u[ijk] + deltaX * u[ijk + 1])
		+ deltaY * ((1 - deltaX) * u[ijk + length + 1] + deltaX * u[ijk + length + 2]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * u[ijk + elemsPerPlane + width] + deltaX * u[ijk + elemsPerPlane + width + 1])
		+ deltaY * ((1 - deltaX) * u[ijk + elemsPerPlane + width + length + 1] + deltaX * u[ijk + elemsPerPlane + width + length + 2]));
}

inline float MACGrid::getV(float x, float y, float z) {
	if (x < 0.5) x = 0.5;
	else if (x > length - 0.5) x = length - 0.5f;
	if (y < 0) y = 0;
	else if (y > width) y = (float)width;
	if (z < 0.5) z = 0.5;
	else if (z > height - 0.5) z = height - 0.5f;
	int i = (int)(x - 0.5f);
	int j = (int)y;
	int k = (int)(z - 0.5f);
	deltaX = x - i - 0.5f;
	deltaY = y - j;
	deltaZ = z - k - 0.5f;
	ijk = t1DV(i, j, k);
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * v[ijk] + deltaX * v[ijk + 1])
		+ deltaY * ((1 - deltaX) * v[ijk + length] + deltaX * v[ijk + length + 1]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * v[ijk + elemsPerPlane + length] + deltaX * v[ijk + elemsPerPlane + length + 1])
		+ deltaY * ((1 - deltaX) * v[ijk + elemsPerPlane + length + length] + deltaX * v[ijk + elemsPerPlane + length + length + 1]));
}

inline float MACGrid::getW(float x, float y, float z) {
	if (x < 0.5) x = 0.5;
	else if (x > length-0.5) x = length - 0.5f;
	if (y < 0.5) y = 0.5;
	else if (y > width - 0.5) y = width - 0.5f;
	if (z < 0) z = 0;
	else if (z > height) z = (float)height;
	int i = (int)(x - 0.5f);
	int j = (int)(y - 0.5f);
	int k = (int)z;
	deltaX = x - i - 0.5f;
	deltaY = y - j - 0.5f;
	deltaZ = z - k;
	ijk = t1DW(i, j, k);
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * w[ijk] + deltaX * w[ijk + 1])
		+ deltaY * ((1 - deltaX) * w[ijk + length] + deltaX * w[ijk + length + 1]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * w[ijk + elemsPerPlane] + deltaX * w[ijk + elemsPerPlane + 1])
		+ deltaY * ((1 - deltaX) * w[ijk + elemsPerPlane + length] + deltaX * w[ijk + elemsPerPlane + length + 1]));
}

inline float MACGrid::interpolate(float* q, float x, float y, float z) {
	if (x < 0) x = 0;
	else if (x > length - 1) x = length - 1.0f;
	if (y < 0) y = 0;
	else if (y > width - 1) y = width - 1.0f;
	if (z < 0) z = 0;
	else if (z > height - 1) z = height - 1.0f;
	int i = int(x);
	int j = int(y);
	int k = int(z);
	deltaX = x - i;
	deltaY = y - j;
	deltaZ = z - k;
	ijk = t1D0(i, j, k);
	return (1 - deltaZ) * ((1 - deltaY) * ((1 - deltaX) * q[ijk] + deltaX * q[ijk + 1])
		+ deltaY * ((1 - deltaX) * q[ijk + length] + deltaX * q[ijk + length + 1]))
		+ deltaZ * ((1 - deltaY) * ((1 - deltaX) * q[ijk + elemsPerPlane] + deltaX * q[ijk + elemsPerPlane + 1])
		+ deltaY * ((1 - deltaX) * q[ijk + elemsPerPlane + length] + deltaX * q[ijk + elemsPerPlane + length + 1]));
}