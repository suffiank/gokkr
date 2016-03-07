CXX=g++
CXXFLAGS=-O2 -std=c++11

OBJECTS=obj/extractkvp.o obj/loadkvp.o obj/timer.o obj/logger.o obj/ylm.o 

all: $(OBJECTS)
	ar rcs obj/util.a $(OBJECTS)

test:
	$(CXX) $(CXXFLAGS) src/main.cc obj/util.a -o bin/test

clean:
	rm -f bin/* obj/*

lib: obj/util.a

obj/util.a: $(OBJECTS)
	ar rcs obj/util.a $(OBJECTS)

obj/%.o: src/%.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

obj/extractkvp.o: src/extractkvp.h
obj/loadkvp.o: src/loadkvp.h
obj/timing.o: src/timing.h
obj/logger.o: src/logger.h
obj/ylm.o: src/ylm.h
