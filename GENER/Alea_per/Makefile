INC_GL = -I /usr/include/i386_linux_gnu
LIB_GL = -lGL -lGLU -lglut

INC_X = -I /usr/include/X11
LIB_X = -L /usr/lib/X11

CXXFLAGS   =  -Wall -O3 -fopenmp -g #-fast -Wefc++ -O3
CXX        = g++ $(INC_GL) $(INC_X)

PROGS= gener3d_coh

all: $(PROGS)

gener3d_coh : gener3d_coh.o
	$(CXX) $^ -o $@ $(CXXFLAGS) $(LIB_GL) $(LIB_X)  -lrt	
		
.cpp.o:
	$(CXX) $(CXXFLAGS)  -c  $<

clean:
	-rm $(PROGS) *.o
