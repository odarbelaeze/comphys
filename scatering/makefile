scatering: scatering.o
	g++ -std=c++11 -O3 -o scatering scatering.o

scatering.o: scatering.cc numerov.h bessel.h
	g++ -std=c++11 -O3 -c scatering.cc

run: scatering
	./scatering
.PHONY: run
