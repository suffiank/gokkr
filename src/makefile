CXX=g++
CXXFLAGS=-O2
LINK= \
  atom/obj/atom.a crystal/obj/crystal.a \
	green0/obj/green0.a util/obj/util.a \
  -L/usr/local/lib -lgsl -lgslcblas -lm


all:
	cd util; make lib
	cd atom; make lib
	cd green0; make lib
	cd crystal; make lib
	$(CXX) $(CXXFLAGS) main.cc $(LINK) -o ../bin/knoxkkr 
