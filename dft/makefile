dft: dft.o
	g++ -o dft dft.o

dft.o: dft.cc numerov.h bisection.h
	g++ -std=c++11 -c dft.cc

run: dft
	./dft

data: dft
	./dft > data

.PHONY: run data
