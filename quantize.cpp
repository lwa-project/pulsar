#include "Python.h"
#include <cmath>
#include <cstdint>

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


template<uint8_t D>
PyObject *OptimizeDataLevels(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *bzero=NULL, *bscale=NULL, *bdata=NULL, *output;
	PyArrayObject *data=NULL, *zeroF=NULL, *scaleF=NULL, *dataF=NULL;
	long ij, i, j, k, nStand, nSamps, nFFT;
	int nChan = 64;
	
	char const* kwlist[] = {"spectra", "nChan", "bzero", "bscale", "bdata", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "Oi|OOO", const_cast<char **>(kwlist), &spectra, &nChan, &bzero, &bscale, &bdata)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 2, 2);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 2-D float32");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan;
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	if( bzero != NULL && bzero != Py_None ) {
		zeroF = (PyArrayObject*) PyArray_ContiguousFromObject(bzero, NPY_FLOAT32, 2, 2);
		if(zeroF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output bzero array to 2-D float32");
			goto fail;
		}
		if(PyArray_DIM(zeroF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "bzero has an unexpected number of stands");
			goto fail;
		}
		if(PyArray_DIM(zeroF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "bzero has an unexpected number of channels");
			goto fail;
		}
	} else {
		zeroF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
		if(zeroF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			goto fail;
		}
	}
	
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	if( bscale != NULL && bscale != Py_None ) {
		scaleF = (PyArrayObject*) PyArray_ContiguousFromObject(bzero, NPY_FLOAT32, 2, 2);
		if(scaleF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output bscale array to 2-D float32");
			goto fail;
		}
		if(PyArray_DIM(scaleF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "bscale has an unexpected number of stands");
			goto fail;
		}
		if(PyArray_DIM(scaleF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "bscale has an unexpected number of channels");
			goto fail;
		}
	} else {
		scaleF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
		if(scaleF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			goto fail;
		}
	}
	
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nSamps;
	if( bdata != NULL && bdata != Py_None ) {
		dataF = (PyArrayObject*) PyArray_ContiguousFromObject(bdata, NPY_UINT8, 2, 2);
		if(dataF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output bdata array to 2-D uint8");
			goto fail;
		}
		if(PyArray_DIM(dataF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "bdata has an unexpected number of stands");
			goto fail;
		}
		if(PyArray_DIM(dataF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "bdata has an unexpected number of channels");
			goto fail;
		}
	} else {
		dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_UINT8, 0);
		if(dataF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			goto fail;
		}
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Go!
	long secStart;
	float tempMin, tempMax, tempV;
	float *a, *b, *c;
	uint8_t *d;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(zeroF);
	c = (float *) PyArray_DATA(scaleF);
	d = (uint8_t *) PyArray_DATA(dataF);
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(secStart, i, j, k, tempV, tempMin, tempMax)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(ij=0; ij<nStand*nChan; ij++) {
			i = ij / nChan;
			j = ij % nChan;
			
			secStart = nSamps*i + j;
			
			k = 0;
			tempMin = *(a + secStart + nChan*k);
			tempMax = *(a + secStart + nChan*k);
			
			for(k=1; k<nFFT; k++) {
				if( *(a + secStart + nChan*k) < tempMin ) {
					tempMin = *(a + secStart + nChan*k);
				} else if( *(a + secStart + nChan*k) > tempMax ) {
					tempMax = *(a + secStart + nChan*k);
				}
			}
			
			*(b + i*nChan + j) = tempMin;
			*(c + i*nChan + j) = (tempMax - tempMin) / (float) ((1<<D) - 1);
			
			for(k=0; k<nFFT; k++) {
				tempV  = *(a + secStart + nChan*k) - tempMin;
				tempV /= *(c + i*nChan + j);
				tempV  = round(tempV);
				
				*(d + secStart + nChan*k) = ((uint8_t) tempV) & ((1<<D) - 1);
			}
		}
	}
	
	Py_END_ALLOW_THREADS
	
	output = Py_BuildValue("(OOO)", PyArray_Return(zeroF), PyArray_Return(scaleF), PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(zeroF);
	Py_XDECREF(scaleF);
	Py_XDECREF(dataF);
	
	return output;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(zeroF);
	Py_XDECREF(scaleF);
	Py_XDECREF(dataF);
	
	return NULL;
}


PyObject *OptimizeDataLevels8Bit(PyObject *self, PyObject *args, PyObject *kwds) {
	return OptimizeDataLevels<8>(self, args, kwds);
}

char OptimizeDataLevels8Bit_doc[] = PyDoc_STR(\
"Given the output of one of the 'Combine' functions, find the bzero and bscale\n\
values that yield the best representation of the data as unsigned characters\n\
and scale the data accordingly.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.float64 (stands by samples) array of combined data\n\
 * LFFT: number of channels per data stream\n\
\n\
Outputs:\n\
 * bzero: 2-D numpy.float32 (stands by channels) array of bzero values\n\
 * bscale: 2-D numpy.float32 (stand by channels) array of bscale values\n\
 * spectra: 2-D numpy.unit8 (stands by samples) array of spectra\n\
");


PyObject *OptimizeDataLevels4Bit(PyObject *self, PyObject *args, PyObject *kwds) {
	return OptimizeDataLevels<4>(self, args, kwds);
}

char OptimizeDataLevels4Bit_doc[] = PyDoc_STR(\
"Given the output of one of the 'Combine' functions, find the bzero and bscale\n\
values that yield the best representation of the data as 4-bit unsigned integers\n\
and scale the data accordingly.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.float64 (stands by samples) array of combined data\n\
 * LFFT: number of channels per data stream\n\
\n\
Outputs:\n\
 * bzero: 2-D numpy.float32 (stands by channels) array of bzero values\n\
 * bscale: 2-D numpy.float32 (stand by channels) array of bscale values\n\
 * spectra: 2-D numpy.unit8 (stands by samples) array of spectra\n\
");
