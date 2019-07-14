# G O K K R  Makefile
#   by Suffian Khan

# compiler, flags, and external linkage
CXX=g++
CXXFLAGS=-O2 -std=c++11 -fopenmp

INCLUDE=-isystem./3rdparty/eigen -isystem./3rdparty/gsl/include
EXT_LINK= -L./3rdparty/gsl/lib -lgsl -lgslcblas -lm

# made use of bash syntax
SHELL=/bin/bash

# source files organized by modules
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

# primary executables and scripts
all: dir obj/util.a obj/atom.a obj/green0.a obj/crystal.a
	@echo " .. linking executable .. "
	$(CXX) $(CXXFLAGS) $(INCLUDE) src/driver.cc $(LINK) -o bin/gokkr

dir:
	@mkdir -p obj/; mkdir -p bin/

clean:
	rm -f bin/* obj/* .depend

print-%:
	@echo '$*=$($*)'

.phony: all clean dir print-%

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
	@if [ -z "$(MSG)" ]; then echo " .. compiling source .. "; fi
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@
	$(eval MSG=DONE)

# external dependencies
3rdparty/gsl:
	@echo " .. downloading, building, and installing 3rdparty/gsl .. "
	@$(eval GSL_BUILD_DIR := $(shell mktemp -d))
	@echo $(GSL_BUILD_DIR)
	@-(curl -o $(GSL_BUILD_DIR)/gsl.tgz http://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz)
	@-(cd $(GSL_BUILD_DIR) && mkdir gsl && tar -xf gsl.tgz -C gsl --strip-components 1)
	@-(cd $(GSL_BUILD_DIR)/gsl && ./configure --prefix=$(CURDIR)/3rdparty/gsl && make && make install)
	@rm -fr $(GSL_BUILD_DIR)

3rdparty/eigen:
	@git submodule init

# object to header dependencies
depend: .depend

.depend: $(SOURCES) 3rdparty/gsl 3rdparty/eigen
	@echo " .. building dependencies .. "
	@rm -f ./.depend
	$(CXX) $(CXXFLAGS) $(INCLUDE) -MM $(SOURCES) >  ./.depend;
	@sed -i 's|\(^.*:\)|obj/\1|' .depend

include .depend
