CXX      = clang++
CXXFLAGS = -Wall -O3 -std=c++11 
CXXFLAGS += -I/home/guterding/local/eigen3/

OBJECTS = main.o files.o random.o solver.o
LDFLAGS = -lboost_system -lboost_filesystem
DEFINES = 

all : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(DEFINES) $(OBJECTS) $(LDFLAGS) -o ctaux
	
main.o : main.cpp solver.hpp random.hpp files.hpp typedefs.hpp 
	$(CXX) $(CXXFLAGS) $(DEFINES) -c main.cpp -o main.o
	
files.o : files.cpp files.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c files.cpp -o files.o
	
random.o : random.cpp random.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c random.cpp -o random.o
	
solver.o : solver.cpp solver.hpp random.hpp files.hpp typedefs.hpp
	$(CXX) $(CXXFLAGS) $(DEFINES) -c solver.cpp -o solver.o

clean : 
	rm ctaux
	rm *.o