CXX=g++
CXXFLAGS=-std=c++11 -fopenmp
INCLUDE=-I/usr/local/Cellar/gsl/1.16/include
LINK=obj/atom.a /usr/local/Cellar/gsl/1.16/lib/libgsl.a

OBJECTS= obj/atom.o

all: bin/atom

atom: 
	$(CXX) $(INCLUDE) $(CXXFLAGS) src/atomscf.cc $(LINK) -o atom
	@mv atom bin/

test:
	$(CXX) $(INCLUDE) $(CXXFLAGS) src/test.cc $(LINK) -o test
	@mv test bin/

clean:
	rm -r bin/* obj/*

lib: obj/atom.a

obj/atom.a: $(OBJECTS)
	$(CXX) -c $(INCLUDE) $(CXXFLAGS) src/atom.cc -o obj/atom.o
	@ar cr obj/atom.a obj/atom.o

obj/%.o: src/%.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

obj/atom.cc: src/atom.h src/src/eqsolver.h
obj/atomscf.cc: src/atom.h ../util/src/mathfn.h
obj/test.cc: src/atom.h ../util/src/mathfn.h
