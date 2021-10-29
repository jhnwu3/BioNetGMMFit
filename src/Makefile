
# CXX to simplify compilation from makefile
CXX = g++

# target linear and nonlinear.
all: snapshot

# target dependencies for main PSO program - also equally messy code inside!
snapshot: main.o fileIO.o linear.o nonlinear.o calc.o
	g++ main.o fileIO.o linear.o nonlinear.o calc.o -o snapshot -fopenmp
main.o: main.cpp main.hpp fileIO.hpp linear.hpp nonlinear.hpp calc.hpp
	g++ -c -O3 -fopenmp main.cpp 
fileIO.o: fileIO.cpp main.hpp fileIO.hpp
	g++ -c fileIO.cpp 
linear.o: linear.cpp main.hpp linear.hpp
	g++ -c linear.cpp -fopenmp
nonlinear.o: nonlinear.cpp main.hpp nonlinear.hpp
	g++ -c nonlinear.cpp -fopenmp
calc.o: calc.cpp main.hpp calc.hpp
	g++ -c calc.cpp

# this target deletes all files produced from the Makefile
# so that a completely new compile of all items is required
clean:
	rm -rf *.o 
