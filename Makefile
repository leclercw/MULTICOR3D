CXXFLAGS   = -Wall -g -O3 -fopenmp #-w -fast -Wefc++ 
CXX        = g++ -w 

SRCS = $(wildcard SOURCE/*/*.cpp)
BIN  = $(wildcard SOURCE/*/*.o)
OBJS= $(SRCS:.cpp=.o)
EXEC= multi3D

all: $(EXEC)

$(EXEC) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^  -lrt	

%.o: %.cpp 
	$(CXX) $(CXXFLAGS) -o $@ -c  $<
	  
SOURCE/MAIN/main.o: $(SRCS)

clean:
	-rm $(BIN)		
	-rm $(EXEC)
