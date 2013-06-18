scatering: scatering.o
	g++ -std=c++11 -o scatering scatering.o

scatering.o: scatering.cc numerov.h bessel.h
	g++ -std=c++11 -c scatering.cc

run: scatering
	./scatering
.PHONY: run
