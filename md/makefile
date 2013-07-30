argon: argon.o
	g++ -std=c++11 -o argon argon.cc

argon.o: argon.h argon.cc
	g++ -std=c++11 -c argon.cc

data: argon
	./argon > data

run: argon
	./argon
.PHONY: run
