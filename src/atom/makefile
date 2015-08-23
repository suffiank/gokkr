CXX=g++
CXXFLAGS=-std=c++11 -fopenmp
INCLUDE=-I/usr/local/Cellar/gsl/1.16/include
LINK=../util/obj/util.a obj/atom.a /usr/local/Cellar/gsl/1.16/lib/libgsl.a

OBJECTS= obj/atom.o

atom: lib
	$(CXX) $(INCLUDE) $(CXXFLAGS) src/atomscf.cc $(LINK) -o atom
	@mv atom bin/

test: lib
	$(CXX) $(INCLUDE) $(CXXFLAGS) src/test_phase.cc $(LINK) -o test_phase
	$(CXX) $(INCLUDE) $(CXXFLAGS) src/test_levels.cc $(LINK) -o test_levels
	@mv test_phase bin/; mv test_levels bin/ 

clean:
	rm -f bin/* obj/*

lib: obj/atom.a

obj/atom.a: $(OBJECTS)
	@ar cr obj/atom.a obj/atom.o

obj/%.o: src/%.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

obj/atom.cc: src/atom.h src/eqsolver.h ../util/src/extractkvp.h
obj/atomscf.cc: src/atom.h ../util/src/mathfn.h ../util/src/loadkvp.h
src/eqsolver.h: ../util/src/mathfn.h
