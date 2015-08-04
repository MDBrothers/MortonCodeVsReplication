#include "plot3d.hpp"

Display                 *dpy;
Window                  root, win;
GLint                   att[]   = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
XVisualInfo             *vi;
GLXContext              glc;
Colormap                cmap;
XSetWindowAttributes    swa;
XWindowAttributes       wa;
XEvent                  xev;
bool terminate = false;

// Type declarations

struct Vertex{
	float coords[4];
	float color[4];
};

struct Matrix4x4{
	float entries[16];
};

// Static globals

static const Matrix4x4 I_MAT_4x4={{1.f,0.f,0.f,0.f,  0.f,1.f,0.f,0.f,  0.f,0.f,1.f,0.f,  0.f,0.f,0.f,1.f}};	  
				
enum buffer {MODEL_VERTICES};
enum object {MODEL};

static Matrix4x4 rotMat = I_MAT_4x4;
static Matrix4x4 transMat = I_MAT_4x4;
static Matrix4x4 projMat = I_MAT_4x4;
float  rot[3]={0.f,0.f,0.f};
float  trans[3]={0.f,0.f,0.f};
float  cam[3]={0.f,0.f,0.f};
float  eye[3]={0.f,0.f,0.f};
float  up[3]={0.f,1.f,0.f};

static unsigned int
   programId,
   vertexShaderId,
   fragmentShaderId,
   modelViewMatLoc,
   projMatLoc,
   buffer[1],
   vao[1];

/* A simple function that will read a file into an allocated char pointer buffer */
char* filetobuf(const std::string file)
{
    FILE *fptr;
    long length;
    char *buf;
 
    fptr = fopen(file.c_str(), "rb"); /* Open file for reading */
    if (!fptr) /* Return NULL on failure */
        return NULL;
    fseek(fptr, 0, SEEK_END); /* Seek to the end of the file */
    length = ftell(fptr); /* Find out how many bytes into the file we are */
    buf = (char*)malloc(length+1); /* Allocate a buffer for the entire length of the file and a null terminator */
    fseek(fptr, 0, SEEK_SET); /* Go back to the beginning of the file */
    fread(buf, length, 1, fptr); /* Read the contents of the file in to the buffer */
    fclose(fptr); /* Close the file */
    buf[length] = 0; /* Null terminator */
 
    return buf; /* Return the buffer */
}

void ExposeFunc(const int numVertices) {
    float  aspect_ratio;
    char   info_string[256];

//Resize the viewport
    XGetWindowAttributes(dpy, win, &wa);
    glViewport(0, 0, wa.width, wa.height);
    aspect_ratio = (float)(wa.width) / (float)(wa.height);

//Update the view and transformation matrices
    for(int i(0); i<3; ++i){
    	rotMat.entries[5*i] = rot[i];
    	//transMat.entries[5*i] = trans[i];
    	//projMat.entries[5*i] = 1.f;
    }

//Draw the model
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawArrays(GL_LINE_STRIP, 0, numVertices);

//Update graphics system stats
    sprintf(info_string, "<up,down,left,right> rotate model * <F1> stop rotation ");

//Swap display buffer with draw buffer
    glXSwapBuffers(dpy, win);
}

//Create a GL window
void CreateWindow() {
    if((dpy = XOpenDisplay(NULL)) == NULL) {
        printf("\n\tcannot connect to x server\n\n");
        exit(0);
    }

    root = DefaultRootWindow(dpy);

    if((vi = glXChooseVisual(dpy, 0, att)) == NULL) {
        printf("\n\tno matching visual\n\n");
        exit(0);
    }

    if((cmap = XCreateColormap(dpy, root, vi->visual, AllocNone)) == 0) {
        printf("\n\tcannot create colormap\n\n");
        exit(0);
    }

    swa.event_mask = KeyPressMask;
    swa.colormap   = cmap;
    win = XCreateWindow(dpy, root, 0, 0, 700, 700, 0, vi->depth, InputOutput, vi->visual, CWColormap | CWEventMask, &swa);
    XStoreName(dpy, win, "OpenGL Animation");
    XMapWindow(dpy, win);
}

//Setup GL context
void SetupGL() {
    char           font_string[128];
    XFontStruct    *font_struct;

//Create GL context and set as current
    if((glc = glXCreateContext(dpy, vi, NULL, GL_TRUE)) == NULL) {
        printf("\n\tcannot create gl context\n\n");
        exit(0);
    }

    glXMakeCurrent(dpy, win, glc);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.f, 0.f, 0.f, 1.f);

    //Locate a font for text display purposes
    for(int font_size = 14; font_size < 32; font_size += 2) {
        sprintf(font_string, "-adobe-courier-*-r-normal--%i-*", font_size);
        font_struct = XLoadQueryFont(dpy, font_string);

        if(font_struct != NULL) {
            glXUseXFont(font_struct->fid, 32, 192, 32);
            break;
        }
    }
}

//Cleanup method
void ExitPlotWindow() {
    glXMakeCurrent(dpy, None, NULL);
    glXDestroyContext(dpy, glc);
    XDestroyWindow(dpy, win);
    XCloseDisplay(dpy);
    exit(0);
}

//User input polling method
void CheckKeyboard() {

    if(XCheckWindowEvent(dpy, win, KeyPressMask, &xev)) {
        char    *key_string = XKeysymToString(XkbKeycodeToKeysym(dpy, xev.xkey.keycode, 0, 0));

        if(strncmp(key_string, "Left", 4) == 0) {
            rot[0] = -45.f;
        }

        else if(strncmp(key_string, "Right", 5) == 0) {
            rot[0] = 45.f;
        }
        else if(strncmp(key_string, "Up", 2) == 0) {
            rot[0] = -45.f;
        }

        else if(strncmp(key_string, "Down", 4) == 0) {
            rot[0] = 45.f;
        }
        else if(strncmp(key_string, "F1", 2) == 0) {
            rot[0] = 0.f;
            rot[0] = 0.f;
        }

        else if(strncmp(key_string, "Escape", 5) == 0) {
            ExitPlotWindow();
	    terminate = true;
        }
    }
}

void plotCurve(const std::valarray<double>& vertices){
	CreateWindow();
	SetupGL();

	// Copy whatever we get into an array of floats
	std::valarray<Vertex> model(vertices.size());
	for(int i = 0; i < vertices.size(); ++i){
		model[i].coords[0] = vertices[i*3];	
		model[i].coords[1] = vertices[i*3+1];	
		model[i].coords[2] = vertices[i*3+2];	
		model[i].coords[3] = 1.f;	
		model[i].color[0] = vertices[i*3];
		model[i].color[1] = vertices[i*3+1];
		model[i].color[2] = vertices[i*3+2];
		model[i].color[3] = 1.f;
	}

	// Compile vertex shader
	char* vertexShader = filetobuf("vertexShader.glsl");
	vertexShaderId = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertexShaderId, 1, &vertexShader, NULL);
	glCompileShader(vertexShaderId);

	// Compile fragment shader
	char* fragmentShader = filetobuf("fragmentShader.glsl");
	fragmentShaderId = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragmentShaderId, 1, &fragmentShader, NULL);
	glCompileShader(fragmentShaderId);

	// Don't need the source code anymore after compilation
	free(fragmentShader);
	free(vertexShader);

	// Put together the programmable part of the graphics pipline, the shaders
	programId = glCreateProgram();
	glAttachShader(programId, vertexShaderId);
	glAttachShader(programId, fragmentShaderId);
	glLinkProgram(programId);
	glUseProgram(programId);

	// Create the Vertex Aray Object and Vertex Buffer Object and associate the data with the vertex shader
	glGenVertexArrays(1, vao);
	glGenBuffers(1, buffer);
	glBindVertexArray(vao[MODEL]);
	glBindBuffer(GL_ARRAY_BUFFER, buffer[MODEL_VERTICES]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex)*model.size(), &model[0], GL_STATIC_DRAW);

	glVertexAttribPointer(0,4,GL_FLOAT,GL_FALSE, sizeof(Vertex), 0);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(1,4,GL_FLOAT,GL_FALSE,sizeof(Vertex), (GLvoid*)sizeof(model[0].coords));
	glEnableVertexAttribArray(1);

	// Upload uniform matrices  
	projMatLoc = glGetUniformLocation(programId,"projMat");
	glUniformMatrix4fv(projMatLoc, 1, GL_TRUE, projMat.entries);

	Matrix4x4 modelViewMat = I_MAT_4x4;	
	modelViewMatLoc = glGetUniformLocation(programId, "modelViewMat");
	glUniformMatrix4fv(modelViewMatLoc, 1, GL_TRUE, modelViewMat.entries);

	while(not terminate){
		ExposeFunc(model.size());
		usleep(1000);
		CheckKeyboard();
	}
}

