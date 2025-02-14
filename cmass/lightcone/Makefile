CC = gcc
CXX = g++

FLAGS = -Wall -Wextra -Wno-reorder -Wno-sign-compare -O3 -ggdb -fPIC
CFLAGS = --std=gnu99 $(FLAGS)
CXXFLAGS = --std=c++17 -fopenmp $(FLAGS)

LFLAGS = -lm -lgsl -lgslcblas -fopenmp

pymangle_src = ./pymangle/pymangle
pymangle = libcmangle.a
healpixlite_src = ./healpix_lite
healpixlite = libhealpixlite.a
cuboidremap_src = cuboidremap/c++
cuboidremap = libcuboidremap.a
lco = lc.o

all: lc


$(lco): lightcone.cpp
	$(CXX) -c $(shell python3 -m pybind11 --includes) -I$(pymangle_src) -I$(healpixlite_src) -I$(cuboidremap_src) -DHAVE_INLINE $(CXXFLAGS) -o $@ $<

lc: $(lco) $(pymangle) $(healpixlite) $(cuboidremap)
	$(CXX) -shared -o lc$(shell python3-config --extension-suffix) $(lco) $(pymangle) $(healpixlite) $(cuboidremap) $(LFLAGS) $(shell python3-config --ldflags) -Wl,-undefined,dynamic_lookup

$(cuboidremap): $(patsubst $(cuboidremap_src)/%.cpp,$(cuboidremap_src)/%.o,$(wildcard $(cuboidremap_src)/*.cpp))
	ar -r $@ $^

$(pymangle): $(patsubst $(pymangle_src)/%.c,$(pymangle_src)/%.o,$(wildcard $(pymangle_src)/*.c))
	ar -r $@ $^

$(healpixlite): $(patsubst $(healpixlite_src)/%.cc,$(healpixlite_src)/%.o,$(wildcard $(healpixlite_src)/*.cc))
	ar -r $@ $^

$(cuboidremap)/%.o: $(cuboidremap_src)/%.cpp
	$(CXX) -c -I$(cuboidremap_src) $(CXXFLAGS) -o $@ $<

$(pymangle_src)/%.o: $(pymangle_src)/%.c
	$(CC) -c -I$(pymangle_src) $(CFLAGS) -o $@ $<

$(healpixlite_src)/%.o: $(healpixlite_src)/%.cc
	$(CXX) -c -I$(healpixlite_src) $(CXXFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -f $(cuboidremap_src)/*.o
	rm -f $(pymangle_src)/*.o
	rm -f $(healpixlite_src)/*.o
	rm -f *.a
	rm -f *.o
	rm -f *.so
