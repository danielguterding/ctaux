CXX      = mpic++.openmpi
CXXFLAGS = -Wall -O3 -I/home/guterding/local/eigen3/

OBJECTS = main.o files.o
LDFLAGS = -lboost_system -lboost_filesystem
DEFINES = 

all : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(DEFINES) $(OBJECTS) $(LDFLAGS) -o ctaux
	
main.o : main.cpp files.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c main.cpp -o main.o
	
files.o : files.cpp files.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c files.cpp -o files.o
	
clean : 
	rm ctaux
	rm *.o