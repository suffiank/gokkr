# G O K K R  Makefile
#   by Suffian Khan

# compiler, flags, and external linkage
# -----------------------------------------------------------------------------
CC=gcc
CXX=g++
CXXFLAGS=-O2 -std=c++11 -fopenmp

INCLUDE= \
  -isystem./3rdparty/eigen \
  -isystem./3rdparty/gsl/include \
  -I./3rdparty/googletest/googletest/include

EXT_LINK= \
  -L./3rdparty/gsl/lib -lgsl -lgslcblas -lm \
  -L./3rdparty/googletest/build/googlemock/gtest -lgtest -lgtest_main

# bake in runpath relative to bin location
EXT_LINK+= \
  -Wl,--enable-new-dtags,-rpath,'$$ORIGIN/../3rdparty/gsl/lib' \
  -Wl,--enable-new-dtags,-rpath,'$$ORIGIN/../3rdparty/googletest/build/lib'

# made use of bash syntax
SHELL=/bin/bash

# check for presence of required executables
# -----------------------------------------------------------------------------
# curl is used to download gnu scientific library
# cmake is ued to build googletest
EXECUTABLES = $(CC) $(CXX) curl cmake
K := $(foreach exec, $(EXECUTABLES),\
   $(if $(shell which $(exec)), dummy, $(error Missing required path executable '$(exec)')))

# source files organized by modules mapped to static library
# -----------------------------------------------------------------------------
UTIL_SRC= extractkvp.cc loadkvp.cc timer.cc logger.cc ylm.cc mathfn.cc
ATOM_SRC= atom.cc
STRC_SRC= export.cc g0block.cc g0aij.cc g0config.cc g0factor.cc \
  g0gaunt.cc g0latt.cc g0matrix.cc faddeeva.cc 
CRYS_SRC= crystal.cc symmetry.cc specialk.cc

# prepend directories to files source, object files
UTIL_SRCS= $(addprefix src/,$(UTIL_SRC))
ATOM_SRCS= $(addprefix src/,$(ATOM_SRC))
STRC_SRCS= $(addprefix src/,$(STRC_SRC))
CRYS_SRCS= $(addprefix src/,$(CRYS_SRC))
SOURCES= $(UTIL_SRCS) $(ATOM_SRCS) $(STRC_SRCS) $(CRYS_SRCS)

UTIL_OBJS= $(addprefix obj/,$(UTIL_SRC:.cc=.o))
ATOM_OBJS= $(addprefix obj/,$(ATOM_SRC:.cc=.o))
STRC_OBJS= $(addprefix obj/,$(STRC_SRC:.cc=.o))
CRYS_OBJS= $(addprefix obj/,$(CRYS_SRC:.cc=.o))
OBJECTS= $(UTIL_OBJS) $(ATOM_OBJS) $(STRC_OBJS) $(CRYS_OBJS)

# add internal static libraries to linking
INT_LINK= $(addprefix obj/,crystal.a atom.a green0.a util.a)
LINK= $(INT_LINK) $(EXT_LINK)

# targets for primary executables, static libraries, and scripts
# -----------------------------------------------------------------------------
all: obj/util.a obj/atom.a obj/green0.a obj/crystal.a src/driver.cc | bin
	@echo " .. linking executable .. "
	$(CXX) $(CXXFLAGS) $(INCLUDE) src/driver.cc $(LINK) -o bin/gokkr

.phony: all clean depend test force print-% 

bin:
	@mkdir -p bin/

test: obj/util.a obj/atom.a obj/green0.a obj/crystal.a src/test_levels.cc force
	@echo " .. building and executing tests .."
	@-(mkdir -p test)
	$(CXX) $(CXXFLAGS) $(INCLUDE) src/test_levels.cc $(LINK) -o test/atomic_solver
	@test/atomic_solver

clean:
	rm -f bin/* obj/* test/* .depend

print-%:
	@echo '$*=$($*)'

force:

# static library dependencies
obj: 
	@mkdir -p obj

obj/util.a: $(UTIL_OBJS)
	ar rcs obj/util.a $(UTIL_OBJS)

obj/crystal.a: $(CRYS_OBJS) 
	ar rcs obj/crystal.a $(CRYS_OBJS)

obj/atom.a: $(ATOM_OBJS) 
	ar rcs obj/atom.a $(ATOM_OBJS)	

obj/green0.a: $(STRC_OBJS)
	ar rcs obj/green0.a $(STRC_OBJS)

# object to source dependencies
obj/%.o: src/%.cc | obj
	@if [ -z "$(MSG)" ]; then echo " .. compiling source .. "; fi
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@
	$(eval MSG=DONE)

# targets for external dependencies
# -----------------------------------------------------------------------------
3rdparty/gsl:
	@echo " .. downloading, building, and installing 3rdparty/gsl .. "
	@$(eval GSL_BUILD_DIR := $(shell mktemp -d))
	@echo $(GSL_BUILD_DIR)
	@-(curl -o $(GSL_BUILD_DIR)/gsl.tgz http://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz)
	@-(cd $(GSL_BUILD_DIR) && mkdir gsl && tar -xf gsl.tgz -C gsl --strip-components 1)
	@-(cd $(GSL_BUILD_DIR)/gsl && ./configure --prefix=$(CURDIR)/3rdparty/gsl && make && make install)
	@rm -fr $(GSL_BUILD_DIR)

3rdparty/eigen:
	@echo " .. initializing 3rdparty/eigen .. "
	@git submodule update --init -- 3rdparty/eigen

3rdparty/googletest/build:
	@echo " .. downloading, building, and installing 3rdparty/googletest .."
	@git submodule update --init -- 3rdparty/googletest
	@mkdir -p 3rdparty/googletest/build && cd 3rdparty/googletest/build && cmake .. && make

# automatic generation of object to header dependencies
# -----------------------------------------------------------------------------
depend: .depend

.depend: $(SOURCES) | 3rdparty/gsl 3rdparty/eigen 3rdparty/googletest/build
	@echo " .. auto-detecting include header dependencies .. "
	@rm -f ./.depend
	$(CXX) $(CXXFLAGS) $(INCLUDE) -MM $(SOURCES) >  ./.depend;
	@sed -i 's|\(^.*:\)|obj/\1|' .depend

-include .depend
