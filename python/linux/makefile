all: _nbody.so

_nbody.so: nbody_wrap.o
	c++ -shared nbody_wrap.o -o _nbody.so -L /lib/python2.7 -lpython2.7

nbody_wrap.o: nbody_wrap.cxx 
	c++ -c nbody_wrap.cxx -o nbody_wrap.o -I /usr/include/python2.7  -L /lib/python2.7

nbody_wrap.cxx: ../nbody.i 
	swig -c++ -python -o ./nbody_wrap.cxx  ../nbody.i

clean:
	rm -rf nbody_wrap.cxx
