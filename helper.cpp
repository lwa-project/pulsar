#include "Python.h"
#include <cmath>

#ifdef _OPENMP
	#include <omp.h>
	
	// OpenMP scheduling method
	#ifndef OMP_SCHEDULER
	#define OMP_SCHEDULER dynamic
	#endif
#endif

#include "numpy/arrayobject.h"


static PyObject *FastAxis0MinMax(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *pulses, *pulsesF;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long i, j, nSamps, nParams;
	
	char const* kwlist[] = {"pulses", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O", const_cast<char **>(kwlist), &pulses)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(pulses, NPY_FLOAT32, 2, 2);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input pulses array to 2-D float32");
		return NULL;
	}
	
	// Get the properties of the data
	nSamps = (long) PyArray_DIM(data, 0);
	nParams = (long) PyArray_DIM(data, 1);
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nParams;
	dims[1] = (npy_intp) 2;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Pointers
	float *a, *b;
	float tempMin, tempMax;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataF);
	
	Py_BEGIN_ALLOW_THREADS
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, tempMin, tempMax)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(j=0; j<nParams; j++) {
			tempMin = 1e200;
			tempMax = -tempMin;
			
			for(i=0; i<nSamps; i++) {
				if( (float) *(a + nParams*i + j) < tempMin ) {
					tempMin = (float) *(a + nParams*i + j);
				} else if( (float) *(a + nParams*i + j) > tempMax ) {
					tempMax = (float) *(a + nParams*i + j);
				}
			}
			
			*(b + 2*j + 0) = (float) tempMin;
			*(b + 2*j + 1) = (float) tempMax;
		}
	}
	
	Py_END_ALLOW_THREADS
	
	Py_XDECREF(data);
	
	pulsesF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);
	
	return pulsesF;
}

PyDoc_STRVAR(FastAxis0MinMax_doc, \
"Given a 2-D numpy.float32 array, compute the minimum and maximum values along\n\
the zeros axis\n\
\n\
Input arguments are:\n\
 * pulses: 2-D numpy.float32 (pulse by parameter) array of data\n\
\n\
Input keywords are:\n\
 None\n\
\n\
Outputs:\n\
 * minmax: 2-D numpy.float32 (parameter by min/max) of the pulses\n\
");


static PyObject *FastHistogram(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *pulses, *edges, *pulsesF;
	PyArrayObject *data=NULL, *bins=NULL, *dataF=NULL;
	
	long i, j, k, l, nSamps, nBins;
	
	char const* kwlist[] = {"pulses", "bins", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO", const_cast<char **>(kwlist), &pulses, &edges)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(pulses, NPY_FLOAT32, 1, 1);
	bins = (PyArrayObject *) PyArray_ContiguousFromObject(edges, NPY_FLOAT32, 1, 1);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input pulses array to 1-D float32");
		return NULL;
	}
	if( bins == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input edges array to 1-D float32");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Get the properties of the data
	nSamps = (long) PyArray_DIM(data, 0);
	
	// Get the properties of the bins
	nBins = (long) PyArray_DIM(bins, 0) - 1;
	
	// Create the output array
	npy_intp dims[1];
	dims[0] = (npy_intp) nBins;
	dataF = (PyArrayObject*) PyArray_ZEROS(1, dims, NPY_INT32, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		Py_XDECREF(bins);
		return NULL;
	}
	
	// Pointers
	float *a, *b;
	int *c;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(bins);
	c = (int *) PyArray_DATA(dataF);
	
	Py_BEGIN_ALLOW_THREADS
	
	// Temporary storage
	int *tempHist;
	for(l=0; l<nBins; l++) {
		*(c + l) = 0;
	}
	
	// Divide up the work
	long nThreads, nSampsPerThread;
	#ifdef _OPENMP
		nThreads = omp_get_max_threads();
		nSampsPerThread = nSamps / nThreads;
		if( nSampsPerThread % nThreads != 0 ) {
			nSampsPerThread += 1;
		}
	#else
		nThreads = 1;
		nSampsPerThread = nSamps;
	#endif
	
// 	printf("Running with %li threads and %li samples/thread\n", nThreads, nSampsPerThread);
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, l, tempHist)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(i=0; i<nThreads; i++) {
			tempHist = (int *) malloc(nBins*sizeof(int));
			for(l=0; l<nBins; l++) {
				*(tempHist + l) = 0;
			}
			
			for(j=0; j<nSampsPerThread; j++) {
				k = i*nSampsPerThread + j;
				if( k >= nSamps ) {
					continue;
				}
				
				for(l=0; l<nBins; l++) {
					if( *(a+k) >= *(b+l) && *(a+k) <= *(b+l+1) ) {
						*(tempHist + l) += 1;
						break;
					}
				}
			}
			
			#ifdef _OPENMP
				#pragma omp critical
			#endif
			{
				for(l=0; l<nBins; l++) {
					*(c + l) += *(tempHist + l);
				}
			}
			
			free(tempHist);
			
		}
	}
	
	Py_END_ALLOW_THREADS
	
	Py_XDECREF(data);
	Py_XDECREF(bins);
	
	pulsesF = Py_BuildValue("OO", PyArray_Return(dataF), edges);
	Py_XDECREF(dataF);
	
	return pulsesF;
}

PyDoc_STRVAR(FastHistogram_doc, \
"Given a 1-D numpy.float32 array of pulses and a 1-D numpy.float32 of\n\
histogram bin eges, compute the histogram.\n\
\n\
Input arguments are:\n\
 * data: 1-D numpy.float32 array of data\n\
 * bins: 1-D numpy.float32 array of bin edges\n\
\n\
Input keywords are:\n\
 None\n\
\n\
Outputs:\n\
 * histogram: 1-D numpy.int32 histogram\n\
");


static PyObject *FastAxis0Mean(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long i, j, k, jk, nStand, nSamps, nChans, iPrime;
	
	if(!PyArg_ParseTuple(args, "O", &spectra)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 3-D float32");
		return NULL;
	}
	
	// Get the properties of the data
	nSamps = (long) PyArray_DIM(data, 0);
	nStand = (long) PyArray_DIM(data, 1);
	nChans = (long) PyArray_DIM(data, 2);
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChans;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Pointers
	float *a, *b;
	float tempV;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataF);
	
	Py_BEGIN_ALLOW_THREADS
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, jk, tempV, iPrime)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(jk=0; jk<nStand*nChans; jk++) {
			j = jk / nChans;
			k = jk % nChans;
			
			tempV = 0.0;
			
			iPrime = 0;
			for(i=0; i<nSamps; i++) {
				if( *(a + nStand*nChans*i + nChans*j + k) == *(a + nStand*nChans*i + nChans*j + k) ) {
					tempV += (float) *(a + nStand*nChans*i + nChans*j + k);
					iPrime++;
				}
			}
			
			if( iPrime > 0 ) {
				tempV /= (float) iPrime;
			} else {
				tempV = 1.0;
			}
			*(b + nChans*j + k) = tempV;
		}
	}
	
	Py_END_ALLOW_THREADS
	
	Py_XDECREF(data);
	
	spectraF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);
	
	return spectraF;
}

PyDoc_STRVAR(FastAxis0Mean_doc, \
"Given a 3-D numpy.float32 array, compute the mean along the zeroth axis\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (time by stands by channels) array of data\n\
\n\
Input keywords are:\n\
 None\n\
\n\
Outputs:\n\
 * mean: 2-D numpy.float32 (stands by channels) of time-averaged spectra\n\
");


static PyObject *FastAxis1MinMax(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long i, j, k, nStand, nSamps, nChans;
	long int chanMin, chanMax;
	chanMin = 0;
	chanMax = -1;
	
	char const* kwlist[] = {"spectra", "chanMin", "chanMax", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|ll", const_cast<char **>(kwlist), &spectra, &chanMin, &chanMax)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 3-D float32");
		return NULL;
	}
	
	// Get the properties of the data
	nSamps = (long) PyArray_DIM(data, 0);
	nStand = (long) PyArray_DIM(data, 1);
	nChans = (long) PyArray_DIM(data, 2);
	if( chanMax < chanMin ) {
		chanMax = nChans;
	}
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) 2;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Pointers
	float *a, *b;
	float tempMin, tempMax;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataF);
	
	Py_BEGIN_ALLOW_THREADS
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, tempMin, tempMax)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(j=0; j<nStand; j++) {
			tempMin = 1e200;
			tempMax = -tempMin;
			
			for(k=chanMin; k<chanMax; k++) {
				for(i=0; i<nSamps; i++) {
					if( (float) *(a + nStand*nChans*i + nChans*j + k) < tempMin ) {
						tempMin = (float) *(a + nStand*nChans*i + nChans*j + k);
					} else if( (float) *(a + nStand*nChans*i + nChans*j + k) > tempMax ) {
						tempMax = (float) *(a + nStand*nChans*i + nChans*j + k);
					}
				}
			}
			
			*(b + 2*j + 0) = (float) tempMin;
			*(b + 2*j + 1) = (float) tempMax;
		}
	}
	
	Py_END_ALLOW_THREADS
	
	Py_XDECREF(data);
	
	spectraF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);
	
	return spectraF;
}

PyDoc_STRVAR(FastAxis1MinMax_doc, \
"Given a 3-D numpy.float32 array, compute the minimum and maximum values along\n\
the first axis\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (time by stands by channels) array of data\n\
\n\
Input keywords are:\n\
 * chanMin: Channel to start at (default = 0)\n\
 * chanMax: Channel to stop at (default = -1 => maximum channel number\n\
\n\
Outputs:\n\
 * minmax: 2-D numpy.float32 (stands by min/max) of the spectra\n\
");


static PyObject *FastAxis0Bandpass(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *bandpass, *spectraF;
	PyArrayObject *data=NULL, *dataB=NULL;
	
	long i, j, k, jk, nStand, nSamps, nChans;
	
	if(!PyArg_ParseTuple(args, "OO", &spectra, &bandpass)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	dataB = (PyArrayObject *) PyArray_ContiguousFromObject(bandpass, NPY_FLOAT32, 2, 2);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 3-D float32");
		return NULL;
	}
	if( dataB == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input bandpass array to 2-D float32");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Get the properties of the data
	nSamps = (long) PyArray_DIM(data, 0);
	nStand = (long) PyArray_DIM(data, 1);
	nChans = (long) PyArray_DIM(data, 2);
	
	// Pointers
	float *a, *b;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataB);
	
	Py_BEGIN_ALLOW_THREADS
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(jk=0; jk<nStand*nChans; jk++) {
			j = jk / nChans;
			k = jk % nChans;
			
			for(i=0; i<nSamps; i++) {
				*(a + nStand*nChans*i + nChans*j + k) /= *(b + nChans*j + k);
			}
		}
	}
	
	Py_END_ALLOW_THREADS
	
	Py_XDECREF(data);
	Py_XDECREF(dataB);
	
	spectraF = Py_BuildValue("i", 1);
	return spectraF;
}

PyDoc_STRVAR(FastAxis0Bandpass_doc, \
"Given a 3-D numpy.float32 array of spectra and a 2-D numpy.float32 of\n\
bandpass shapes, apply a bandpass correction to the spectra.\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (time by stands by channels) array of data\n\
 * bandpass: 2-D numpy.float32 (stands by channels) array of bandpass shapes\n\
\n\
Input keywords are:\n\
 None\n\
\n\
Outputs:\n\
 * bandpassed: 3-D numpy.float32 (time by stands by channels) of bandpass-\n\
               corrected spectra\n\
");


int cmpfloat(const void *a, const void *b) {
	/*
	 * cmpfloat - Comparison function for qsort-ing an array of float values.
	 */
	if( *(const float *) a < *(const float *) b ) {
		return -1;
	}
	return *(const float *) a > *(const float *) b;
}


static PyObject *FastAxis0Median(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data=NULL, *dataF=NULL;
	
	long i, j, k, nStand, nSamps, nChans, iPrime;
	
	if(!PyArg_ParseTuple(args, "O", &spectra)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 3-D float32");
		return NULL;
	}
	
	// Get the properties of the data
	nSamps = (long) PyArray_DIM(data, 0);
	nStand = (long) PyArray_DIM(data, 1);
	nChans = (long) PyArray_DIM(data, 2);
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[2];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChans;
	dataF = (PyArrayObject*) PyArray_ZEROS(2, dims, NPY_FLOAT32, 0);
	if(dataF == NULL) {
		PyErr_Format(PyExc_MemoryError, "Cannot create output array");
		Py_XDECREF(data);
		return NULL;
	}
	
	// Pointers
	float *a, *b;
	a = (float *) PyArray_DATA(data);
	b = (float *) PyArray_DATA(dataF);
	
	float *tempV;
	
	Py_BEGIN_ALLOW_THREADS
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, tempV, iPrime)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(k=0; k<nChans; k++) {
			tempV = (float *) malloc(nSamps*sizeof(float));
			
			for(j=0; j<nStand; j++) {
				
				iPrime = 0;
				for(i=0; i<nSamps; i++) {
					*(tempV + iPrime) = *(a + nStand*nChans*i + nChans*j + k);
					if( *(tempV + iPrime) == *(tempV + iPrime) ) {
						iPrime++;
					}
				}
				
				if( iPrime > 0 ) {
					qsort(tempV, iPrime, sizeof(float), cmpfloat);
					
					*(b + nChans*j + k) = *(tempV + iPrime/2);
				} else {
					*(b + nChans*j + k) = 1.0;
				}
			}
			
			free(tempV);
		}
	}
	
	Py_END_ALLOW_THREADS
	
	Py_XDECREF(data);
	
	spectraF = Py_BuildValue("O", PyArray_Return(dataF));
	Py_XDECREF(dataF);
	
	return spectraF;
}

PyDoc_STRVAR(FastAxis0Median_doc, \
"Given a 3-D numpy.float32 array, compute the approximate median along the\n\
zeroth axis\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (time by stands by channels) array of data\n\
\n\
Input keywords are:\n\
 None\n\
\n\
Outputs:\n\
 * median: 2-D numpy.float32 (stands by channels) of the median spectra\n\
\n\
.. note::\n\
\tThis function uses a median-of-medians method and may give slightly\n\
\tdifferent results relative to numpy.median() for arrays with even\n\
\tnumbers of elements.\n\
");


static PyObject *FastAxis1Percentiles5And99(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *spectra, *spectraF;
	PyArrayObject *data=NULL;
	
	long i, j, k, nStand, nSamps, nChans, iPrime;
	long int stand, chanMin, chanMax;
	stand = 0;
	chanMin = 0;
	chanMax = -1;
	
	char const* kwlist[] = {"spectra", "stand", "chanMin", "chanMax", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "Ol|ll", const_cast<char **>(kwlist), &spectra, &stand, &chanMin, &chanMax)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(spectra, NPY_FLOAT32, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectra array to 3-D float32");
		return NULL;
	}
	
	// Get the properties of the data
	nSamps = (long) PyArray_DIM(data, 0);
	nStand = (long) PyArray_DIM(data, 1);
	nChans = (long) PyArray_DIM(data, 2);
	if( chanMax < chanMin ) {
		chanMax = nChans;
	}
	
	// Pointers
	float *temp5, *temp99;
	
	Py_BEGIN_ALLOW_THREADS
	
	// Pointers
	float *a, *tempV;
	a = (float *) PyArray_DATA(data);
	
	temp5  = (float *) malloc((chanMax-chanMin)*sizeof(float));
	temp99 = (float *) malloc((chanMax-chanMin)*sizeof(float));
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(i, j, k, tempV, iPrime)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(k=chanMin; k<chanMax; k++) {
			j = stand;
			
			tempV = (float *) malloc(nSamps*sizeof(float));
			
			iPrime = 0;
			for(i=0; i<nSamps; i++) {
				*(tempV + iPrime) = *(a + nStand*nChans*i + nChans*j + k);
				if( *(tempV + iPrime) == *(tempV + iPrime) ){
					iPrime++;
				}
			}
			
			qsort(tempV, iPrime, sizeof(float), cmpfloat);
			
			if( iPrime > 0 ) {
				*(temp5  + k - chanMin) = *(tempV + (int) (iPrime * 0.05));
				*(temp99 + k - chanMin) = *(tempV + (int) (iPrime * 0.99));
			} else {
				*(temp5  + k - chanMin) = 0.1;
				*(temp99 + k - chanMin) = 0.2;
			}
			
			free(tempV);
		}
	}
	
	qsort(temp5,  (chanMax-chanMin), sizeof(float), cmpfloat);
	qsort(temp99, (chanMax-chanMin), sizeof(float), cmpfloat);
	
	Py_END_ALLOW_THREADS
	
	Py_XDECREF(data);
	
	spectraF = Py_BuildValue("ff", *(temp5 + (int) ((chanMax-chanMin)*0.05)), *(temp99 + (int) ((chanMax-chanMin) * 0.99)));
	
	free(temp5);
	free(temp99);
	
	return spectraF;
}

PyDoc_STRVAR(FastAxis1Percentiles5And99_doc, \
"Given a 3-D numpy.float32 array, compute the approximate fifth and 99th \n\
percentiles for a particular index along the first axis.\n\
\n\
Input arguments are:\n\
 * spectra: 3-D numpy.float32 (time by stands by channels) array of data\n\
 * stand: index along the first axis to compute these values\n\
\n\
Input keywords are:\n\
 * chanMin: Channel to start at (default = 0)\n\
 * chanMax: Channel to stop at (default = -1 => maximum channel number\n\
\n\
Outputs:\n\
 * percentiles: two-element tuple of the fifth and 99th percentile values\n\
\n\
.. note::\n\
\tThis function uses percentile-of-percentile method to compute\n\
\tthe returned values.  These should be *close* to those returned\n\
\tby the scipy.stats.scoreatpercentile function for many cases.\n\
");


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef HelperMethods[] = {
	{"FastAxis0MinMax",            (PyCFunction) FastAxis0MinMax,            METH_VARARGS|METH_KEYWORDS, FastAxis0MinMax_doc}, 
	{"FastHistogram",              (PyCFunction) FastHistogram,              METH_VARARGS|METH_KEYWORDS, FastHistogram_doc}, 
	{"FastAxis0Mean",              (PyCFunction) FastAxis0Mean,              METH_VARARGS,               FastAxis0Mean_doc}, 
	{"FastAxis1MinMax",            (PyCFunction) FastAxis1MinMax,            METH_VARARGS|METH_KEYWORDS, FastAxis1MinMax_doc}, 
	{"FastAxis0Bandpass",          (PyCFunction) FastAxis0Bandpass,          METH_VARARGS,               FastAxis0Bandpass_doc}, 
	{"FastAxis0Median",            (PyCFunction) FastAxis0Median,            METH_VARARGS,               FastAxis0Median_doc}, 
	{"FastAxis1Percentiles5And99", (PyCFunction) FastAxis1Percentiles5And99, METH_VARARGS|METH_KEYWORDS, FastAxis1Percentiles5And99_doc},
	{NULL,                         NULL,                                     0,                          NULL}
};

PyDoc_STRVAR(helper_doc, \
"HELPER - The Heuristic Extension aLlowing Parallel Efficiency in Rendering\n\
\n\
Extension to help plotHDF.py along by running the data computations in\n\
parallel.  See the individual functions for more details.");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC PyInit__helper(void) {
	PyObject *m, *all;

	// Module definitions and functions
	static struct PyModuleDef moduledef = {
		 PyModuleDef_HEAD_INIT, "_helper", helper_doc, -1, HelperMethods
	};
	m  = PyModule_Create(&moduledef);
	if( m == NULL ) {
        return NULL;
    }
	import_array();
	
	// Version information
	PyModule_AddObject(m, "__version__", PyUnicode_FromString("0.1"));
	
	// Function listings
	all = PyList_New(0);
	PyList_Append(all, PyUnicode_FromString("FastAxis0MinMax"));
	PyList_Append(all, PyUnicode_FromString("FastHistogram"));
	PyList_Append(all, PyUnicode_FromString("FastAxis0Mean"));
	PyList_Append(all, PyUnicode_FromString("FastAxis1MinMax"));
	PyList_Append(all, PyUnicode_FromString("FastAxis0Bandpass"));
	PyList_Append(all, PyUnicode_FromString("FastAxis0Median"));
	PyList_Append(all, PyUnicode_FromString("FastAxis1Percentiles5And99"));
	PyModule_AddObject(m, "__all__", all);
	
	return m;
}
