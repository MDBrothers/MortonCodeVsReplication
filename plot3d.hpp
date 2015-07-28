#ifndef PLOT3D_H
#define PLOT3D_H

#include<X11/Xlib.h>
#include<X11/XKBlib.h>

#include<GL/gl.h>
#include<GL/glx.h>
#include<GL/glext.h>
#include<GL/glu.h>

#include<valarray>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
#include<iostream>
#include<sstream>
#include<string>
#include<exception>
#include<algorithm>
#include <cstdlib>

void plotCurve(const std::valarray<double>& vertices); 

#endif
