CXX=mpic++

#CXXFLAGS= -std=c++17 -Wall -pedantic -march=native -ffast-math -fopenmp
#CXXFLAGS= -std=c++17 -O3 -ffast-math -funroll-all-loops -fprefetch-loop-arrays -fopenmp
CXXFLAGS= -std=c++17 -O3 -fopenmp -funsafe-math-optimizations -floop-parallelize-all -fprefetch-loop-arrays -funroll-all-loops -ffast-math -fexpensive-optimizations

.DEFAULT_GOAL := all

all: parallel

parallel: Makefile parallel.cpp

	$(CXX) parallel.cpp -o parallel -v $(CXXFLAGS)


.PHONY: clean
clean:
	rm parallel
