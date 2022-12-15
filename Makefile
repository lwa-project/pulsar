PYTHON ?= python

.PHONY: all
all: _psr.so _helper.so

_psr.so: psr.cpp utils.cpp fft.cpp kurtosis.cpp dedispersion.cpp reduce.cpp quantize.cpp setup.py
	$(PYTHON) setup.py build
	mv build/lib*/*.so .
	rm -rf build
	
_helper.so: helper.cpp setup.py
	$(PYTHON) setup.py build
	mv build/lib*/*.so .
	rm -rf build
	
clean:
	rm -rf _psr.*so _psr.*dylib _helper.*so _helper.*dylib
