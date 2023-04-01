#include "Python.h"
#include <cmath>
#include <complex>
#include <fftw3.h>

#ifdef _OPENMP
	#include <omp.h>
	
	// OpenMP scheduling method
	#ifndef OMP_SCHEDULER
	#define OMP_SCHEDULER dynamic
	#endif
#endif

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL psr_ARRAY_API
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

#include "psr.hpp"


PyObject *ComputeSKMask(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF;
	PyArrayObject *data=NULL, *dataF=NULL;
	double lower, upper;
	long ij, i, j, k, nStand, nSamps, nChan, nFFT;
	
	if(!PyArg_ParseTuple(args, "Odd", &signals, &lower, &upper)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 3-D complex64");
		return NULL;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nChan  = (long) PyArray_DIM(data, 1);
	nFFT   = (long) PyArray_DIM(data, 2);
	
	// Find out how large the output array needs to be and initialize it
	nSamps = nChan*nFFT;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Go!
	long secStart;
	float tempV, tempV2, temp2V;
	Complex32 *a;
	float *b;
	a = (Complex32 *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataF);
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(secStart, i, j, k, tempV, tempV2, temp2V)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(ij=0; ij<nStand*nChan; ij++) {
			i = ij / nChan;
			j = ij % nChan;
			
			secStart = nSamps*i + nFFT*j;
			
			tempV2 = 0.0;
			temp2V = 0.0;
			for(k=0; k<nFFT; k++) {
				tempV  = abs2(*(a + secStart + k));
				temp2V += tempV*tempV;
				tempV2 += tempV;
			}
			
			tempV  = nFFT*temp2V / (tempV2*tempV2) - 1.0;
			tempV *= (nFFT + 1.0)/(nFFT - 1.0);
			
			if( tempV < lower || tempV > upper ) {
				*(b + nChan*i + j) = 0.0;
			} else {
				*(b + nChan*i + j) = 1.0;
			}
		}
	}
	
	Py_END_ALLOW_THREADS
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return signalsF;
}

char ComputeSKMask_doc[] = PyDoc_STR(\
"Given the output of PulsarEngineRaw, calculate channel weights based on\n\
spectral kurtosis\n\
\n\
Input arguments are:\n\
 * signals: 3-D numpy.complex64 (stands by channels by integrations) array\n\
   of data\n\
 * lower: lower spectral kurtosis limit\n\
 * upper: upper spectral kurtosis limit\n\
\n\
Outputs:\n\
 * weight: 2-D numpy.float32 (stands by channels) of data weights\n\
");


PyObject *ComputePseudoSKMask(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF;
	PyArrayObject *data=NULL, *dataF=NULL;
	double lower, upper;
	long ij, i, j, k, nStand, nSamps, nChan, nFFT, skN;
	
	if(!PyArg_ParseTuple(args, "Olldd", &signals, &nChan, &skN, &lower, &upper)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_DOUBLE, 2, 2);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 2-D float64");
		return NULL;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps  = (long) PyArray_DIM(data, 1);
	nFFT = nSamps / nChan;
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Go!
	long secStart;
	double tempV, temp2V, tempV2;
	double *a;
	float *b;
	a = (double *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataF);
	
	#ifdef _OPENMP
		omp_set_dynamic(0);
		#pragma omp parallel default(shared) private(secStart, i, j, k, tempV, temp2V, tempV2)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(ij=0; ij<nStand*nChan; ij++) {
			i = ij / nChan;
			j = ij % nChan;
			
			secStart = nSamps*i;
			
			temp2V = 0.0;
			tempV2 = 0.0;
			for(k=0; k<nFFT; k++) {
				tempV  = *(a + secStart + nChan*k + j);

				temp2V += tempV*tempV;
				tempV2 += tempV;
			}
			
			tempV  = nFFT*temp2V / (tempV2*tempV2) - 1.0;
			tempV *= (nFFT*skN + 1.0)/(nFFT - 1.0);
			
			if( tempV < lower || tempV > upper ) {
				*(b + nChan*i + j) = 0.0;
			} else {
				*(b + nChan*i + j) = 1.0;
			}
		}
	}
	
	Py_END_ALLOW_THREADS
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return signalsF;
}

char ComputePseudoSKMask_doc[] = PyDoc_STR(\
"Given the output of DR spectrometer, calculate channel weights based on\n\
psuedo-spectral kurtosis\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.complex64 (stands by channels/integrations) array\n\
   of data to from DR spectrometer\n\
 * LFFT: FFT length\n\
 * N: number of FFT windows per integration\n\
 * lower: lower spectral kurtosis limit\n\
 * upper: upper spectral kurtosis limit\n\
\n\
Outputs:\n\
 * weight: 2-D numpy.float32 (stands by channels) of data weights\n\
");
