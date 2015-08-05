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

static unsigned int
   programId,
   vertexShaderId,
   fragmentShaderId,
   modelViewMatLoc,
   projMatLoc,
   rotXMatLoc,
   rotYMatLoc,
   rotZMatLoc,
   transMatLoc,
   buffer[1],
   vao[1],
   plotStyle(0);


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

static Matrix4x4 rotXMat = I_MAT_4x4;
static Matrix4x4 rotYMat = I_MAT_4x4;
static Matrix4x4 rotZMat = I_MAT_4x4;
static Matrix4x4 transMat = I_MAT_4x4;
static Matrix4x4 projMat = I_MAT_4x4;

float  rot[3]={0.f,0.f,0.f};
float  trans[3]={0.f,0.f,-5.f};
float  cam[3]={0.f,0.f,0.f};
float  eye[3]={0.f,0.f,0.f};
float  up[3]={0.f,1.f,0.f};

void updateXRot(Matrix4x4 &rotMat, float rot){
	rotMat.entries[5] = cos(rot);
	rotMat.entries[6] = -sin(rot);
	rotMat.entries[9] = sin(rot);
	rotMat.entries[10] = cos(rot);
}

void updateYRot(Matrix4x4 &rotMat, float rot){
	rotMat.entries[0] = cos(rot);
	rotMat.entries[2] = sin(rot);
	rotMat.entries[8] = -sin(rot);
	rotMat.entries[10] = cos(rot);
}

void updateZRot(Matrix4x4 &rotMat, float rot){
	rotMat.entries[0] = cos(rot);
	rotMat.entries[1] = -sin(rot);
	rotMat.entries[4] = sin(rot);
	rotMat.entries[5] = cos(rot);
}

void updateTransMat(Matrix4x4 &myTransMat, float *trans){
	myTransMat.entries[3] = trans[0];
	myTransMat.entries[7] = trans[1];
	myTransMat.entries[11] = trans[2];
}

void buildPerspProjMat(Matrix4x4& myProjMat, float fov, float aspect, float znear, float zfar){
	float xymax = znear * tan(fov * M_PI / 360.f);
	float ymin = -xymax;
	float xmin = -xymax;

	float width = xymax - xmin;
	float height = xymax - ymin;

	float depth = zfar - znear;
	float q = -(zfar + znear) / depth;
	float qn = -2.f * (zfar * znear) / depth;

	float w = 2.f * znear / width;
	w = w / aspect;
	float h = 2.f * znear / height;

	myProjMat.entries[0]  = w;
	myProjMat.entries[5]  = h;
	myProjMat.entries[10] = q;
	myProjMat.entries[11] = -1.f;
	myProjMat.entries[14] = qn;
	myProjMat.entries[15] = 0.f;
}

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
    updateXRot(rotXMat, rot[0]);
    updateYRot(rotYMat, rot[1]);
    updateZRot(rotZMat, rot[2]);
    updateTransMat(transMat, trans);
    buildPerspProjMat(projMat, 60.f, aspect_ratio, 1.f, 1000.f);

//Draw the model
    glPointSize(3.0);
    glLineWidth(2.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(plotStyle == 0) glDrawArrays(GL_LINE_STRIP, 0, numVertices);
    else glDrawArrays(GL_POINTS, 0, numVertices);



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

	//std::cout << XKeysymToString(XkbKeycodeToKeysym(dpy, xev.xkey.keycode, 0, 0)) << std::endl;
 

        if(strncmp(key_string, "KP_Home", 7) == 0) {
            rot[0] += 5.f/180.f*M_PI;
        }
        else if(strncmp(key_string, "KP_Begin", 8) == 0) {
            rot[1] += 5.f/180.f*M_PI;
        }
        else if(strncmp(key_string, "KP_Next", 7) == 0) {
            rot[2] += 5.f/180.f*M_PI;
        }
        else if(strncmp(key_string, "w", 1) == 0) {
	    trans[1] += 1.E-1;
	}
	else if(strncmp(key_string, "space", 5) == 0){
	    plotStyle += 1;
	    plotStyle %= 2;
	}
        else if(strncmp(key_string, "s", 1) == 0) {
	    trans[1] -= 1.E-1;
	}
        else if(strncmp(key_string, "a", 1) == 0) {
	    trans[0] -= 1.E-1;
	}
        else if(strncmp(key_string, "d", 1) == 0) {
	    trans[0] += 1.E-1;
	}
	else if(strncmp(key_string, "e", 1) == 0) {
	    trans[2] += 1.E-1;
	}
	else if(strncmp(key_string, "q", 1) == 0) {
	    trans[2] -= 1.E-1;
	}
        else if(strncmp(key_string, "F1", 2) == 0) {
            rot[1] = 0.f;
            rot[2] = 0.f;
            rot[3] = 0.f;
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

	// Get locations for uniform matrices  
	Matrix4x4 modelViewMat = I_MAT_4x4;	
	projMatLoc = glGetUniformLocation(programId,"projMat");
	rotXMatLoc = glGetUniformLocation(programId,"rotXMat");
	rotYMatLoc = glGetUniformLocation(programId,"rotYMat");
	rotZMatLoc = glGetUniformLocation(programId,"rotZMat");
	transMatLoc = glGetUniformLocation(programId, "transMat");
	modelViewMatLoc = glGetUniformLocation(programId, "modelViewMat");

	while(not terminate){
		// Upload update uniform matrices
		glUniformMatrix4fv(projMatLoc, 1, GL_TRUE, projMat.entries);
		glUniformMatrix4fv(rotXMatLoc, 1, GL_TRUE, rotXMat.entries);
		glUniformMatrix4fv(rotYMatLoc, 1, GL_TRUE, rotYMat.entries);
		glUniformMatrix4fv(rotZMatLoc, 1, GL_TRUE, rotZMat.entries);
		glUniformMatrix4fv(transMatLoc, 1, GL_TRUE, transMat.entries);
		glUniformMatrix4fv(modelViewMatLoc, 1, GL_TRUE, modelViewMat.entries);
		// Change window size if needed and draw
		ExposeFunc(model.size());
		usleep(1000);
		CheckKeyboard();
	}
}

