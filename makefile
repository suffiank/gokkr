# gokkr Makefile
#   by Suffian Khan

# compiler, flags, and link-line
CXX=g++
CXXFLAGS=-O2 -std=c++11 -fopenmp

INCLUDE=-isystem../../tools/eigen-3.2.8/
EXT_LINK= -L/usr/local/lib -lgsl -lgslcblas -lm

# all object files sorted by module
UTIL_SRC= extractkvp.cc loadkvp.cc timer.cc logger.cc ylm.cc mathfn.cc
ATOM_SRC= atom.cc
STRC_SRC= export.cc g0block.cc g0aij.cc g0config.cc g0factor.cc \
  g0gaunt.cc g0latt.cc g0matrix.cc faddeeva.cc 
CRYS_SRC= crystal.cc symmetry.cc specialk.cc

UTIL_SRCS= $(addprefix src/,$(UTIL_SRC))
ATOM_SRCS= $(addprefix src/,$(ATOM_SRC))
STRC_SRCS= $(addprefix src/,$(STRC_SRC))
CRYS_SRCS= $(addprefix src/,$(CRYS_SRC))
SOURCES= $(UTIL_SRCS) $(ATOM_SRCS) $(STRC_SRCS) $(CRYS_SRCS)

UTIL_OBJS= $(addprefix obj/,$(UTIL_SRC:.cc=.o))
ATOM_OBJS= $(addprefix obj/,$(ATOM_SRC:.cc=.o))
STRC_OBJS= $(addprefix obj/,$(STRC_SRC:.cc=.o))
CRYS_OBJS= $(addprefix obj/,$(CRYS_SRC:.cc=.o))
OBJECTS= $(UTIL_OBJS) $(ATOM_OBJS) $(STRC_OBJCS) $(CRYS_OBJS)

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

# object to header dependencies
depend: .depend

.depend: $(SOURCES)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^ -MF  ./.depend;

include .depend

# obj/extractkvp.o: src/extractkvp.h
# obj/loadkvp.o: src/loadkvp.h
# obj/timing.o: src/timing.h
# obj/logger.o: src/logger.h
# obj/ylm.o: src/ylm.h
# obj/mathfn.o: src/mathfn.h
# obj/crystal.o: src/crystal.h
# obj/symmetry.o: src/crystal.h
# obj/specialk.o: src/crystal.h
# src/crystal.h: src/symmetry.h src/specialk.h src/mat3.h
# obj/g0*.o: src/green0.h
# obj/export.o: src/green0.h
# obj/g0config.o: src/extractkvp.h
