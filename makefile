all: parse.o
	g++ -Wall -O4 -DDEBUG -static -std=c++11 -pthread -o parse parse.o lib.a
parse.o: parse.cpp header.h
	g++ -Wall -O4 -DDEBUG -static -std=c++11 -o parse.o -c parse.cpp