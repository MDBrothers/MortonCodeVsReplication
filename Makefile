CXX_FLAGS = -std=c++11 -O4
THIRD_PARTY = -I/usr/local/hdf5/include -I/usr/local/trilinos/include -I/usr/local/Trilinos/include -L/usr/local/hdf5/lib -L/usr/local/trilinos/lib -L/usr/local/Trilinos/lib -lflann -lmpi -lhdf5 -lepetra -lepetraext -lteuchoscore -lteuchoscomm -lteuchosparameterlist -lnox -lnoxepetra -lloca -llocaepetra -I/usr/local/boost/include

compile_tests: 
	data_test

libplot3d.so: plot3d.cpp plot3d.hpp vertexShader.glsl fragmentShader.glsl
	mpicxx $(CXX_FLAGS) -shared -fPIC plot3d.cpp -o libplot3d.so -lGLU -lGL -DGL_GLEXT_PROTOTYPES

libdata.so: data.cpp data.hpp
	mpicxx $(CXX_FLAGS) -shared -fPIC data.cpp -o libdata.so $(THIRD_PARTY)

clean:
	rm *.so main *_test
