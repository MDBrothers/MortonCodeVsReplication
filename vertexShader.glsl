#version 430 // Specify which version of GLSL we are using.

layout(location=0) in vec4 modelCoords;
layout(location=1) in vec4 modelColors;

uniform mat4 projMat;
uniform mat4 modelViewMat;
uniform mat4 rotXMat, rotYMat, rotZMat, transMat;

out vec4 colorsExport;

void main(void) 
{ 
    gl_Position = projMat * modelViewMat * transMat * rotXMat * rotYMat * rotZMat * modelCoords;
    colorsExport = modelColors + vec4(.1,.1,.1,0.f);
}
