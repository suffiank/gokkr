# gokkr Makefile
#   by Suffian Khan

# compiler, flags, and link-line
CXX=g++
CXXFLAGS=-O2 -std=c++11 -fopenmp

INCLUDE=-isystem../../util/eigen-3.2.8/
EXT_LINK= -L/usr/local/lib -lgsl -lgslcblas -lm

# all object files sorted by module
UTIL_OBJ= extractkvp.o loadkvp.o timer.o logger.o ylm.o mathfn.o
ATOM_OBJ= atom.o
STRC_OBJ= export.o g0block.o g0aij.o g0config.o g0factor.o \
  g0gaunt.o g0latt.o g0matrix.o faddeeva.o 
CRYS_OBJ= crystal.o symmetry.o specialk.o

UTIL_OBJS= $(addprefix obj/,$(UTIL_OBJ))
ATOM_OBJS= $(addprefix obj/,$(ATOM_OBJ))
STRC_OBJS= $(addprefix obj/,$(STRC_OBJ))
CRYS_OBJS= $(addprefix obj/,$(CRYS_OBJ))

INT_LINK= $(addprefix obj/,util.a atom.a green0.a crystal.a)
LINK= $(INT_LINK) $(EXT_LINK)

# primary executables and clean script
all: direct obj/util.a obj/atom.a obj/green0.a obj/crystal.a
	$(CXX) $(CXXFLAGS) $(INCLUDE) src/driver.cc $(LINK) -o bin/gokkr

direct:
	mkdir -p obj/; mkdir -p bin/

clean:
	rm -f bin/* obj/*

.phony: all clean direct

# static library dependencies
obj/util.a: $(UTIL_OBJS)
	ar rcs obj/util.a $(UTIL_OBJS)

obj/crystal.a: $(CRYS_OBJS)
	ar rcs obj/crystal.a $(CRYS_OBJS)

obj/atom.a: $(ATOM_OBJS)
	ar rcs obj/atom.a $(ATOM_OBJS)	

obj/green0.a: $(STRC_OBJS)
	ar rcs obj/green0.a $(STRC_OBJS)

# object to source dependencies
obj/%.o: src/%.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

obj/%.o: ext/%.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

# object to source dependencies
obj/extractkvp.o: src/extractkvp.h
obj/loadkvp.o: src/loadkvp.h
obj/timing.o: src/timing.h
obj/logger.o: src/logger.h
obj/ylm.o: src/ylm.h
obj/mathfn.o: src/mathfn.h
obj/crystal.o: src/crystal.h
obj/symmetry.o: src/crystal.h
obj/specialk.o: src/crystal.h
src/crystal.h: src/symmetry.h src/specialk.h src/mat3.h
obj/g0*.o: src/green0.h
obj/export.o: src/green0.h
obj/g0config.o: src/extractkvp.h
