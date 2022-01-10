
# CXX to simplify compilation from makefile
CXX = g++

# target linear and nonlinear.
all: sig

# target dependencies for main PSO program - also equally messy code inside!
sig: fileIO.o linear.o nonlinear.o calc.o main.o system.o
	g++ main.o fileIO.o linear.o nonlinear.o calc.o system.o -o sig -fopenmp
calc.o: calc.cpp main.hpp calc.hpp
	g++ -c calc.cpp
fileIO.o: fileIO.cpp main.hpp fileIO.hpp 
	g++ -c fileIO.cpp 
linear.o: linear.cpp main.hpp linear.hpp calc.hpp
	g++ -c linear.cpp -fopenmp
nonlinear.o: nonlinear.cpp main.hpp nonlinear.hpp calc.hpp
	g++ -c nonlinear.cpp -fopenmp
system.o: ../system.cpp system.hpp main.hpp
	g++ -c ../system.cpp 
main.o: main.cpp main.hpp fileIO.hpp linear.hpp nonlinear.hpp system.hpp calc.hpp
	g++ -c -O3 -fopenmp main.cpp 


# this target deletes all files produced from the Makefile
# so that a completely new compile of all items is required
clean:
	rm -rf *.o 
