#version 430 // Specify which version of GLSL we are using.

in vec4 colorsExport;

out vec4 colorsOut;

void main(void) 
{
	colorsOut = colorsExport;
}
