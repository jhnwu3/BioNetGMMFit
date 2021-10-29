
# CXX to simplify compilation from makefile
CXX = g++

# target linear and nonlinear.
all: linear nonlinear input

linear: linear.o
	g++ linear.o -o linear -fopenmp
linear.o: linear.cpp
	g++ -c -O3 linear.cpp -o linear.o -fopenmp

	
nonlinear: nonlinear.o
	g++ nonlinear.o -o nonlinear -fopenmp
nonlinear.o: nonlinear.cpp
	g++ -c -O3 nonlinear.cpp -o nonlinear.o -fopenmp

input: input.o
	g++ input.o -o input -fopenmp
input.o: input.cpp
	g++ -c -O3 input.cpp -o input.o -fopenmp


# this target deletes all files produced from the Makefile
# so that a completely new compile of all items is required
clean:
	rm -rf *.o 
