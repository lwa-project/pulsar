PYTHON ?= python

.PHONY: all
all: _psr.so _helper.so drx2drxi

_psr.so: psr.cpp utils.cpp fft.cpp kurtosis.cpp dedispersion.cpp reduce.cpp quantize.cpp setup.py
	$(PYTHON) setup.py build
	mv build/lib*/*.so .
	rm -rf build
	
_helper.so: helper.cpp setup.py
	$(PYTHON) setup.py build
	mv build/lib*/*.so .
	rm -rf build

drx.o: drx.cpp drx.hpp lwa.hpp
	g++ -c -o drx.o drx.cpp

drxi.o: drxi.cpp drxi.hpp drx.hpp
	g++ -c -o drxi.o drxi.cpp

drx2drxi.o: drx2drxi.cpp
	g++ -c -o drx2drxi.o drx2drxi.cpp

drx2drxi: drx2drxi.o drx.o drxi.o
	g++ -o drx2drxi drx2drxi.o drx.o drxi.o

clean:
	rm -rf _psr.*so _psr.*dylib _helper.*so _helper.*dylib drx2drxi drx.o drxi.o drx2drxi.o
