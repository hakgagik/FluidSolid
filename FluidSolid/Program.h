#ifndef _GLUT_
#define _GLUT_
#include "Dependencies/glew/glew.h"
#include "Dependencies/freeglut/freeglut.h"
#endif

#include <string>
#include <unordered_map>

class Program {
public:
	static std::unordered_map<std::string, GLuint> programs;

	static std::string readShader(const char*);
	static GLuint createProgram(const char*, const char*, std::string);
	static GLuint getProgram(std::string);
};