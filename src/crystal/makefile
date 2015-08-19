CXX=g++
CXXFLAGS=-O2
LINK=obj/crystal.a

OBJECTS= \
	obj/crystal.o obj/symmetry.o obj/specialk.o

test: obj/crystal.a
	$(CXX) $(CXXFLAGS) src/main.cc $(LINK) -o bin/test 

clean:
	rm -r bin/* obj/*

lib: obj/crystal.a

obj/crystal.a: $(OBJECTS)
	ar rcs obj/crystal.a $(OBJECTS)

obj/%.o: src/%.cc
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

obj/crystal.o: src/crystal.h
obj/symmetry.o: src/crystal.h
obj/specialk.o: src/crystal.h
src/crystal.h: src/symmetry.h src/specialk.h ../util/src/mat3.h
