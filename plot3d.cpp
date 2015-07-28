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
bool terminate = false, on_gpu = false;

float                   TimeCounter, LastFrameTimeCounter, DT, prevTime = 0.0, FPS;
struct timeval          tv, tv0;
int                     Frame = 1, FramesPerFPS;

GLfloat                 rotation_matrix[16];
float                   rot_z = 0.0, rot_y = 0.0;

/* A simple function that will read a file into an allocated char pointer buffer */
char* filetobuf(char *file)
{
    FILE *fptr;
    long length;
    char *buf;
 
    fptr = fopen(file, "rb"); /* Open file for reading */
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

void RotateModel() {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glRotatef(.01, 0.0, 1.0, 0.0);
    glRotatef(.01, 0.0, 0.0, 1.0);
    glMultMatrixf(rotation_matrix);
    glGetFloatv(GL_MODELVIEW_MATRIX, rotation_matrix);
}

void ExposeFunc(const int numVertices) {
    float  aspect_ratio;
    char   info_string[256];

//Resize the viewport
    XGetWindowAttributes(dpy, win, &wa);
    glViewport(0, 0, wa.width, wa.height);
    aspect_ratio = (float)(wa.width) / (float)(wa.height);

//Set up projection (orthographic vs. perspective) and model view (camera location and aim)
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
// glOrtho(-2.50*aspect_ratio, 2.50*aspect_ratio, -2.50, 2.50, -200.0, 200.0);
    gluPerspective(50.0, aspect_ratio, 1.0, 400.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    float eyeX(float(0.0)), eyeY(float(0.0)), eyeZ(-300.0);
    float aimX(eyeX), aimY(eyeY), aimZ(0.0);
    float upX(0.f), upY(1.f), upZ(0.f);

    gluLookAt(eyeX, eyeY, eyeZ, aimX, aimY, aimZ, upX, upY, upZ);
    glMultMatrixf(rotation_matrix);

//Draw the model
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawArrays(GL_LINE_STRIP, 0, numVertices);

//Display graphics system stats
/*
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, (float)wa.width, 0, (float)wa.height, -1., 1.);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glColor3f(1.0, 1.0, 1.0);

    sprintf(info_string, "%4.1f seconds * %4.1f fps at %i x %i", TimeCounter, FPS, wa.width, wa.height);
    glRasterPos2i(10, 10);
    glCallLists(strlen(info_string), GL_UNSIGNED_BYTE, info_string);

    sprintf(info_string, "<up,down,left,right> rotate model * <F1> stop rotation ");
    glRasterPos2i(10, wa.height-32);
    glCallLists(strlen(info_string), GL_UNSIGNED_BYTE, info_string);
*/

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
    glClearColor(0.00, 0.00, 0.00, 1.00);

    //Locate a font for text display purposes
    for(int font_size = 14; font_size < 32; font_size += 2) {
        sprintf(font_string, "-adobe-courier-*-r-normal--%i-*", font_size);
        font_struct = XLoadQueryFont(dpy, font_string);

        if(font_struct != NULL) {
            glXUseXFont(font_struct->fid, 32, 192, 32);
            break;
        }
    }

//Initialize the rotation matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glGetFloatv(GL_MODELVIEW_MATRIX, rotation_matrix);
}

//Time counter functions
void InitTimeCounter() {
    gettimeofday(&tv0, NULL);
    FramesPerFPS = 30;
}

void UpdateTimeCounter() {
    LastFrameTimeCounter = TimeCounter;
    gettimeofday(&tv, NULL);
    TimeCounter = (float)(tv.tv_sec-tv0.tv_sec) + 0.000001*((float)(tv.tv_usec-tv0.tv_usec));
    DT = TimeCounter - LastFrameTimeCounter;
}

void CalculateFPS() {
    Frame ++;

    if((Frame%FramesPerFPS) == 0) {
        FPS = ((float)(FramesPerFPS)) / (TimeCounter-prevTime);
        prevTime = TimeCounter;
    }
}

//Cleanup method
void ExitPlotWindow() {
    glXMakeCurrent(dpy, None, NULL);
    glXDestroyContext(dpy, glc);
    XDestroyWindow(dpy, win);
    XCloseDisplay(dpy);
    //exit(0);
}

//User input polling method
void CheckKeyboard() {

    if(XCheckWindowEvent(dpy, win, KeyPressMask, &xev)) {
        char    *key_string = XKeysymToString(XkbKeycodeToKeysym(dpy, xev.xkey.keycode, 0, 0));

        if(strncmp(key_string, "Left", 4) == 0) {
            rot_z = -45.0;
        }

        else if(strncmp(key_string, "Right", 5) == 0) {
            rot_z = 45.0;
        }
        else if(strncmp(key_string, "Up", 2) == 0) {
            rot_y = -45.0;
        }

        else if(strncmp(key_string, "Down", 4) == 0) {
            rot_y = 45.0;
        }
        else if(strncmp(key_string, "F1", 2) == 0) {
            rot_y = 0.0;
            rot_z = 0.0;
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
	InitTimeCounter();

	// Copy whatever we get into an array of floats
	std::valarray<float> verticesFormatted(vertices.size());
	for(int i = 0; i < vertices.size(); ++i){
		verticesFormatted[i] = vertices[i];	
	}		
	const int numVertices = verticesFormatted.size()/3;

	/*
 	* Copy the vertex values into a buffer on the GPU for fast drawing
 	*/

	/* Create a variable to hold the VBO identifier */
	GLuint curveVBO;
	 
	/* This is a handle to the shader program */
	GLuint shaderProgram;

	/* These pointers will receive the contents of our shader source code files */
	GLchar *vertexSource, *fragmentSource;

	/* These are handles used to reference the shaders */
	GLuint vertexShader, fragmentShader;

	const unsigned int shaderAttribute = 0;

	/*---------------------- Initialise VBO - (Note: do only once, at start of program) ---------------------*/
	/* Create a new VBO and use the variable "curveVBO" to store the VBO id */
	glGenBuffers(1, &curveVBO);

	/* Make the new VBO active */
	glBindBuffer(GL_ARRAY_BUFFER, curveVBO);

	/* Upload vertex data to the video device */
	glBufferData(GL_ARRAY_BUFFER, verticesFormatted.size()*sizeof(float), &verticesFormatted[0], GL_STATIC_DRAW);

	/* Specify that our coordinate data is going into attribute index 0(shaderAttribute), and contains three floats per vertex */
	glVertexAttribPointer(shaderAttribute, 3, GL_FLOAT, GL_FALSE, 0, 0);

	/* Enable attribute index 0(shaderAttribute) as being used */
	glEnableVertexAttribArray(shaderAttribute);

	/* Make the new VBO active. */
	glBindBuffer(GL_ARRAY_BUFFER, curveVBO);
	/*-------------------------------------------------------------------------------------------------------*/

	/*--------------------- Load Vertex and Fragment shaders from files and compile them --------------------*/
	/* Read our shaders into the appropriate buffers */
	vertexSource = filetobuf("plot3dVertexShader.vert");
	fragmentSource = filetobuf("plot3dFragmentShader.frag");

	/* Assign our handles a "name" to new shader objects */
	vertexShader = glCreateShader(GL_VERTEX_SHADER);
	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

	/* Associate the source code buffers with each handle */
	glShaderSource(vertexShader, 1, (const GLchar**)&vertexSource, 0);
	glShaderSource(fragmentShader, 1, (const GLchar**)&fragmentSource, 0);

	/* Free the temporary allocated memory */
	free(vertexSource);
	free(fragmentSource);

	/* Compile our shader objects */
	glCompileShader(vertexShader);
	glCompileShader(fragmentShader);
	/*-------------------------------------------------------------------------------------------------------*/

	/*-------------------- Create shader program, attach shaders to it and then link it ---------------------*/
	/* Assign our program handle a "name" */
	shaderProgram = glCreateProgram();

	/* Attach our shaders to our program */
	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);

	/* Bind attribute index 0 (shaderAttribute) to in_Position*/
	/* "in_Position" will represent "verticesFormatted" array's contents in the vertex shader */
	glBindAttribLocation(shaderProgram, shaderAttribute, "in_Position");

	/* Link shader program*/
	glLinkProgram(shaderProgram);
	/*-------------------------------------------------------------------------------------------------------*/

	/* Set shader program as being actively used */
	glUseProgram(shaderProgram);

	while(not terminate){
		UpdateTimeCounter();
		CalculateFPS();
		RotateModel();
		ExposeFunc(numVertices);
		usleep(1000);
		CheckKeyboard();
	}
}

