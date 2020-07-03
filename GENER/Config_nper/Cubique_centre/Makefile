INC_GL = -I /usr/include/i386_linux_gnu
LIB_GL = -lGL -lGLU -lglut

INC_X = -I /usr/include/X11
LIB_X = -L /usr/lib/X11

CXX        = g++ $(INC_GL) $(INC_X)

PROGS= gener

all: $(PROGS)

gener : gener.o
	$(CXX) $^ -o $@ $(LIB_GL) $(LIB_X)  -lrt
	
.cpp.o:
	$(CXX) -c  $<

clean:
	-rm $(PROGS) *.o
