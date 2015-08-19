CXX=g++
CXXFLAGS=-O2 -std=c++11

OBJECTS=obj/extractkvp.o obj/loadkvp.o obj/timing.o 

all: $(OBJECTS)
	ar rcs obj/util.a $(OBJECTS)

clean:
	rm -r bin/* obj/*

lib: obj/util.a

obj/util.a: $(OBJECTS)
	ar rcs obj/util.a $(OBJECTS)

obj/%.o: src/%.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

obj/extractkvp.o: src/extractkvp.h
obj/loadkvp.o: src/loadkvp.h
obj/timing.o: src/timing.h
