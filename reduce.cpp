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


void Linear(long nStand,
            long nSamps,
            Complex32 const *x,
            Complex32 const *y,
            float* out) {
	for(long i=0; i<nStand; i+=2) {
		Complex32 aX, aY;
		aX = *(x + i*nSamps);
		aY = *(y + i*nSamps);
		*(out + i*nSamps) = abs2(aX);		// XX
		*(out + (i+1)*nSamps) = abs2(aY);		// YY
	}
}


void Intensity(long nStand,
               long nSamps,
               Complex32 const *x,
               Complex32 const *y,
               float* out) {
	for(long i=0; i<nStand; i+=2) {
		Complex32 aX, aY;
		aX = *(x + i*nSamps);
		aY = *(y + i*nSamps);
		*(out + (i/2)*nSamps) = abs2(aX) + abs2(aY);	// XX + YY
	}
}


void Circular(long nStand,
              long nSamps,
              Complex32 const *x,
              Complex32 const *y,
              float* out) {
	for(long i=0; i<nStand; i+=2) {
		Complex32 aX, aYI;
		aX = *(x + i*nSamps);
		aYI = *(y + i*nSamps) * Complex32(0,1);
		*(out + i*nSamps) = abs2(aX + aYI) / 2.0;	// LL
		*(out + (i+1)*nSamps) = abs2(aX - aYI) / 2.0;	// RR
	}
}


void Stokes(long nStand,
            long nSamps,
            Complex32 const *x,
            Complex32 const *y,
            float* out) {
	for(long i=0; i<nStand; i+=2) {
		Complex32 aX, aY, aUV;
		aX = *(x + i*nSamps);
		aY = *(y + i*nSamps);
		aUV = aX * conj(aY);
		*(out + ((i/2)*4 + 0)*nSamps) = abs2(aX) + abs2(aY);	// XX + YY = I
		*(out + ((i/2)*4 + 1)*nSamps) = abs2(aX) - abs2(aY);	// XX + YY = Q
		*(out + ((i/2)*4 + 2)*nSamps) = 2.0*real(aUV);	// U
		*(out + ((i/2)*4 + 3)*nSamps) = -2.0*imag(aUV);		// V
	}
}


typedef void (*ReductionOp)(long, long, Complex32 const*, Complex32 const*, float*);


template<ReductionOp R>
void reduce_engine(long nStand,
	                 long nChan,
									 long nFFT,
									 Complex32 const* data,
									 float* reduced) {
  // Setup
	long nSamps = nChan*nFFT;
	int jk, j, k;
	
	Py_BEGIN_ALLOW_THREADS
	
	// Go!
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(j)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(jk=0; jk<nSamps; jk++) {
			j = jk / nFFT;
			k = jk % nFFT;
		  
			R(nStand, nSamps, (data+nFFT*j+k), (data+nSamps+nFFT*j+k), (reduced+nChan*k+j));
		}
	}
	
	Py_END_ALLOW_THREADS
}


PyObject *CombineToIntensity(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF=NULL;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long nStand, nSamps, nChan, nFFT;
	
	char const* kwlist[] = {"signals", "signalsF", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", const_cast<char **>(kwlist), &signals, &signalsF)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 3-D complex64");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nChan  = (long) PyArray_DIM(data, 1);
	nFFT   = (long) PyArray_DIM(data, 2);
	
	// Find out how large the output array needs to be and initialize it
	nSamps = nChan*nFFT;
	npy_intp dims[2];
	dims[0] = (npy_intp) (nStand/2);
	dims[1] = (npy_intp) nSamps;
	if( signalsF != NULL && signalsF != Py_None ) {
		dataF = (PyArrayObject *) PyArray_ContiguousFromObject(signalsF, NPY_FLOAT32, 2, 2);
		if(dataF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output signalsF array to 2-D float32");
			goto fail;
		}
		if(PyArray_DIM(dataF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of stands");
			goto fail;
		}
		if(PyArray_DIM(dataF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of samples");
			goto fail;
		}
	} else {
		dataF = (PyArrayObject*) PyArray_EMPTY(2, dims, NPY_FLOAT32, 0);
		if(dataF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			goto fail;
		}
	}
	
	reduce_engine<Intensity>(nStand, nChan, nFFT,
		                       (Complex32 const*) PyArray_DATA(data),
													 (float*) PyArray_DATA(dataF));
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return signalsF;
	
fail:
  Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return NULL;
}

char CombineToIntensity_doc[] = PyDoc_STR(\
"Given the output of PulsarEngineRaw, calculate total intensity spectra\n\
\n\
Input arguments are:\n\
 * signals: 3-D numpy.complex64 (stands by channels by integrations) array\n\
   of data to FFT\n\
\n\
Outputs:\n\
 * sub-integration: 2-D numpy.float64 (stands by channels) of spectra data\n\
");


PyObject *CombineToLinear(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF=NULL;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long nStand, nSamps, nChan, nFFT;
	
	char const* kwlist[] = {"signals", "signalsF", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", const_cast<char **>(kwlist), &signals, &signalsF)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 3-D complex64");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nChan  = (long) PyArray_DIM(data, 1);
	nFFT   = (long) PyArray_DIM(data, 2);
	
	// Find out how large the output array needs to be and initialize it
	nSamps = nChan*nFFT;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nSamps;
	if( signalsF != NULL && signalsF != Py_None ) {
		dataF = (PyArrayObject *) PyArray_ContiguousFromObject(signalsF, NPY_FLOAT32, 2, 2);
		if(dataF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output signalsF array to 2-D float32");
			goto fail;
		}
		if(PyArray_DIM(dataF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of stands");
			goto fail;
		}
		if(PyArray_DIM(dataF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of samples");
			goto fail;
		}
	} else {
		dataF = (PyArrayObject*) PyArray_EMPTY(2, dims, NPY_FLOAT32, 0);
		if(dataF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			goto fail;
		}
	}
	
	reduce_engine<Linear>(nStand, nChan, nFFT,
		                    (Complex32 const*) PyArray_DATA(data),
												(float*) PyArray_DATA(dataF));
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return signalsF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return NULL;
}

char CombineToLinear_doc[] = PyDoc_STR(\
"Given the output of PulsarEngineRaw, calculate linear polarization spectra\n\
\n\
Input arguments are:\n\
 * signals: 3-D numpy.complex64 (stands by channels by integrations) array\n\
   of data to FFT\n\
\n\
Outputs:\n\
 * sub-integration: 2-D numpy.float64 (stands by channels) of spectra data\n\
");


PyObject *CombineToCircular(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF=NULL;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long nStand, nSamps, nChan, nFFT;
	
	char const* kwlist[] = {"signals", "signalsF", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", const_cast<char **>(kwlist), &signals, &signalsF)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 3-D complex64");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nChan  = (long) PyArray_DIM(data, 1);
	nFFT   = (long) PyArray_DIM(data, 2);
	
	// Find out how large the output array needs to be and initialize it
	nSamps = nChan*nFFT;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nSamps;
	if( signalsF != NULL && signalsF != Py_None ) {
		dataF = (PyArrayObject *) PyArray_ContiguousFromObject(signalsF, NPY_FLOAT32, 2, 2);
		if(dataF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output signalsF array to 2-D float32");
			goto fail;
		}
		if(PyArray_DIM(dataF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of stands");
			goto fail;
		}
		if(PyArray_DIM(dataF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of samples");
			goto fail;
		}
	} else {
		dataF = (PyArrayObject*) PyArray_EMPTY(2, dims, NPY_FLOAT32, 0);
		if(dataF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			goto fail;
		}
	}
	
	reduce_engine<Circular>(nStand, nChan, nFFT,
		                      (Complex32 const*) PyArray_DATA(data),
													(float*) PyArray_DATA(dataF));
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return signalsF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(dataF);

	return NULL;
}

char CombineToCircular_doc[] = PyDoc_STR(\
"Given the output of PulsarEngineRaw, calculate circular polarization spectra\n\
\n\
Input arguments are:\n\
 * signals: 3-D numpy.complex64 (stands by channels by integrations) array\n\
   of data to FFT\n\
\n\
Outputs:\n\
 * sub-integration: 2-D numpy.float64 (stands by channels) of spectra data\n\
");


PyObject *CombineToStokes(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF=NULL;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long nStand, nSamps, nChan, nFFT;
	
	char const* kwlist[] = {"signals", "signalsF", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", const_cast<char **>(kwlist), &signals, &signalsF)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 3-D complex64");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nChan  = (long) PyArray_DIM(data, 1);
	nFFT   = (long) PyArray_DIM(data, 2);
	
	// Find out how large the output array needs to be and initialize it
	nSamps = nChan*nFFT;
	npy_intp dims[2];
	dims[0] = (npy_intp) (4*(nStand/2));
	dims[1] = (npy_intp) nSamps;
	if( signalsF != NULL && signalsF != Py_None ) {
		dataF = (PyArrayObject *) PyArray_ContiguousFromObject(signalsF, NPY_FLOAT32, 2, 2);
		if(dataF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output signalsF array to 2-D float32");
			goto fail;
		}
		if(PyArray_DIM(dataF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of stands");
			goto fail;
		}
		if(PyArray_DIM(dataF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of samples");
			goto fail;
		}
	} else {
		dataF = (PyArrayObject*) PyArray_EMPTY(2, dims, NPY_FLOAT32, 0);
		if(dataF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			goto fail;
		}
	}
	
	reduce_engine<Stokes>(nStand, nChan, nFFT,
		                    (Complex32 const*) PyArray_DATA(data),
											  (float*) PyArray_DATA(dataF));
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(dataF);

	return signalsF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return NULL;
}

char CombineToStokes_doc[] = PyDoc_STR(\
"Given the output of PulsarEngineRaw, calculate Stokes parameter spectra\n\
\n\
Input arguments are:\n\
 * signals: 3-D numpy.complex64 (stands by channels by integrations) array\n\
   of data to FFT\n\
\n\
Outputs:\n\
 * sub-integration: 2-D numpy.float64 (stands by channels) of spectra data\n\
");
