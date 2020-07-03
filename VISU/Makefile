# Makefile d�velopper par fortin pour la visualisation de MULTICOR SOUS OPENGL
# modifi� le 4 mars 2006

SHELL = /bin/sh

CC = g++ -O3 -w
EXEC=multicorvisu
SRCS= $(wildcard */*.cpp)
BIN = $(wildcard */*.o)
OBJS= $(SRCS:%.cpp=%.o)

INCLUDES = #-I/INCL
LIBSDIR = -L/usr/X11R6/lib
LIBS = -lGL -lGLU -lglut  -lX11 #-lXi -lm -lXmu  

# Cibles 
all : $(EXEC)

clean: 
	-rm $(EXEC)
	-rm $(BIN)

$(EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LIBSDIR) $(LIBS)

%.o: %.c
	$(CC) -c -o $@ $< $(INCLUDES)


