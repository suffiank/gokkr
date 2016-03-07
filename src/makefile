CXX=g++
CXXFLAGS=-O2 -std=c++11 -fopenmp
LINK= \
  atom/obj/atom.a crystal/obj/crystal.a \
	green0/obj/green0.a util/obj/util.a \
  -L/usr/local/lib -lgsl -lgslcblas -lm
INCLUDE=-I../../../util/eigen-3.2.5/ \
	-Iatom/src/ -Igreen0/src -Iutil/src \
	-Icrystal/src/

all:
	cd util; make lib
	cd atom; make lib
	cd green0; make lib
	cd crystal; make lib
	$(CXX) $(CXXFLAGS) $(INCLUDE) driver.cc $(LINK) -o ../bin/knoxkkr

clean:
	rm -f green0/obj/*
	rm -f crystal/obj/*
	rm -f util/obj/*
	rm -f atom/obj/*
