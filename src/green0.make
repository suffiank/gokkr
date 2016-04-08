CXX=g++
CXXFLAGS=-O2 -std=c++11
INCLUDE=-I/usr/local/Cellar/gsl/1.16/include
LINK= obj/green0.a ../util/obj/util.a -L/usr/local/lib -lgsl -lgslcblas -lm
EXEC=bin/strconst

OBJECTS= \
	obj/export.o obj/g0block.o obj/g0aij.o \
	obj/g0config.o obj/g0factor.o obj/g0gaunt.o \
	obj/g0latt.o obj/g0matrix.o obj/faddeeva.o

all: obj/green0.a ../util/src/loadkvp.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) src/main.cc $(LINK) -o $(EXEC)

clean:
	rm -f bin/* obj/*

lib: obj/green0.a

obj/green0.a: $(OBJECTS)
	ar rcs obj/green0.a $(OBJECTS)

obj/%.o: src/%.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

obj/%.o: ext/%.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@


