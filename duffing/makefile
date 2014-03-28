CPPFLAGS = -std=c++0x -O3

duff: duff.o duffing.o 
	g++ $(CPPFLAGS) -o duff duffing.o duff.o

duffing.o: duffing.cc duffing.h
	g++ -c $(CPPFLAGS) duffing.cc

duff.o: duff.cc
	g++ -c $(CPPFLAGS) duff.cc
