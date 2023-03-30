#pragma once

#include <complex>

/*
 Complex types
*/

typedef std::complex<float> Complex32;
typedef std::complex<double> Complex64;

/*
 Macro for 2*pi
*/

#define TPI (2*NPY_PI*Complex64(0,1))


/*
 Dispersion constant in MHz^2 s / pc cm^-3
*/

#define DCONST (double) 4.148808e3


/* 
 Support Functions
*/

inline char* PyUnicode_AsString(PyObject *ob) {
    PyObject *enc;
    char *cstr;
    enc = PyUnicode_AsEncodedString(ob, "utf-8", "Error");
    if( enc == NULL ) {
        PyErr_Format(PyExc_ValueError, "Cannot encode string");
        return NULL;
    }
    cstr = PyBytes_AsString(enc);
    Py_XDECREF(enc);
    return cstr;
}

// fft.c
void read_wisdom(char*, PyObject*);
// utils.c
extern PyObject *BindToCore(PyObject*, PyObject*, PyObject*);
extern char BindToCore_doc[];
extern PyObject *BindOpenMPToCores(PyObject*, PyObject*, PyObject*);
extern char BindOpenMPToCores_doc[];

/*
  Complex magnitude squared functions
*/

template<typename T>
inline T abs2(std::complex<T> z) {
	return z.real()*z.real() + z.imag()*z.imag();
}


/*
 FFT Functions
*/

// fft.c
extern PyObject *PulsarEngineRaw(PyObject*, PyObject*, PyObject*);
extern char PulsarEngineRaw_doc[];
extern PyObject *PulsarEngineRawWindow(PyObject*, PyObject*, PyObject*);
extern char PulsarEngineRawWindow_doc[];
extern PyObject *PhaseRotator(PyObject*, PyObject*, PyObject*);
extern char PhaseRotator_doc[];

/*
 Spectral Kurtosis (RFI Flagging) Functions
*/

// kurtosis.c
extern PyObject *ComputeSKMask(PyObject*, PyObject*, PyObject*);
extern char ComputeSKMask_doc[];
extern PyObject *ComputePseudoSKMask(PyObject*, PyObject*, PyObject*);
extern char ComputePseudoSKMask_doc[];


/*
 Coherent Dedispersion
*/

// dedispersion.c
extern PyObject *MultiChannelCD(PyObject*, PyObject*, PyObject*);
extern char MultiChannelCD_doc[];


/*
 Reduction Functions
*/

// reduce.c
extern PyObject *CombineToIntensity(PyObject*, PyObject*, PyObject*);
extern char CombineToIntensity_doc[];
extern PyObject *CombineToLinear(PyObject*, PyObject*, PyObject*);
extern char CombineToLinear_doc[];
extern PyObject *CombineToCircular(PyObject*, PyObject*, PyObject*);
extern char CombineToCircular_doc[];
extern PyObject *CombineToStokes(PyObject*, PyObject*, PyObject*);
extern char CombineToStokes_doc[];


/*
 Quantizing Functions
*/

// quantize.c
extern PyObject *OptimizeDataLevels8Bit(PyObject*, PyObject*, PyObject*);
extern char OptimizeDataLevels8Bit_doc[];
extern PyObject *OptimizeDataLevels4Bit(PyObject*, PyObject*, PyObject*);
extern char OptimizeDataLevels4Bit_doc[];
