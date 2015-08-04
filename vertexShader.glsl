#version 430 // Specify which version of GLSL we are using.

layout(location=0) in vec4 modelCoords;
layout(location=1) in vec4 modelColors;

uniform mat4 projMat;
uniform mat4 modelViewMat;

out vec4 colorsExport;

void main(void) 
{
    gl_Position = projMat * modelViewMat * modelCoords;
    colorsExport = modelColors;
}
