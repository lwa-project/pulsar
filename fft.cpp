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

#include "psr.h"


/*
 Load in FFTW wisdom.  Based on the read_wisdom function in PRESTO.
*/

void read_wisdom(char *filename, PyObject *m) {
	int status = 0;
	FILE *wisdomfile;
	
	wisdomfile = fopen(filename, "r");
	if( wisdomfile != NULL ) {
		status = fftwf_import_wisdom_from_file(wisdomfile);
		PyModule_AddObject(m, "useWisdom", PyBool_FromLong(status));
		fclose(wisdomfile);
	} else {
		PyModule_AddObject(m, "useWisdom", PyBool_FromLong(status));
	}
}


template<typename InType, typename OutType>
void pulsar_engine(long nStand,
									 long nSamps,
									 long nFFT,
									 int nChan,
									 InType const* data,
									 double const* window,
									 OutType* fdomain) {
	// Setup
  long ij, i, j, k;
	
	Py_BEGIN_ALLOW_THREADS
	
	// Create the FFTW plan
	Complex32 *inP, *in;
	inP = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan);
	fftwf_plan p;
	p = fftwf_plan_dft_1d(nChan,
		                    reinterpret_cast<fftwf_complex*>(inP),
												reinterpret_cast<fftwf_complex*>(inP),
												FFTW_FORWARD, FFTW_ESTIMATE);
	
	// FFT
	long secStart;
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(in, secStart, i, j, k)
	#endif
	{
		in = (Complex32*) fftwf_malloc(sizeof(Complex32) * nChan);
		
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(ij=0; ij<nStand*nFFT; ij++) {
			i = ij / nFFT;
			j = ij % nFFT;
			
			secStart = nSamps * i + nChan*j;
			
			for(k=0; k<nChan; k++) {
				in[k]  = Complex32(*(data + 2*secStart + 2*k + 0), *(data + 2*secStart + 2*k + 1));
				if( window != NULL ) {
          in[k] *= *(window + k);
        }
			}
			
			fftwf_execute_dft(p,
											  reinterpret_cast<fftwf_complex*>(in),
												reinterpret_cast<fftwf_complex*>(in));
			
			for(k=0; k<nChan/2+nChan%2; k++) {
				*(fdomain + nFFT*nChan*i + nFFT*(k + nChan/2) + j) = in[k] / (float) sqrt(nChan);
			}
			for(k=nChan/2+nChan%2; k<nChan; k++) {
				*(fdomain + nFFT*nChan*i + nFFT*(k - nChan/2 - nChan%2) + j) = in[k] / (float) sqrt(nChan);
			}
		}
		
		fftwf_free(in);
	}
	fftwf_destroy_plan(p);
	fftwf_free(inP);
	
	Py_END_ALLOW_THREADS
}

PyObject *PulsarEngineRaw(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *signalsF=NULL;
	PyArrayObject *data=NULL, *dataF=NULL;
	int nChan = 64;
	
	long nStand, nSamps, nFFT;
	
	char const* kwlist[] = {"signals", "LFFT", "signalsF", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|iO", const_cast<char **>(kwlist), &signals, &nChan, &signalsF)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 2-D complex64");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dims[2] = (npy_intp) (nSamps/nChan);
	if( signalsF != NULL && signalsF != Py_None ) {
		dataF = (PyArrayObject *) PyArray_ContiguousFromObject(signalsF, NPY_COMPLEX64, 3, 3);
		if(dataF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output signalsF array to 3-D complex64");
		  goto fail;
		}
		if(PyArray_DIM(dataF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of stands");
			goto fail;
		}
		if(PyArray_DIM(dataF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of channels");
			goto fail;
		}
		if(PyArray_DIM(dataF, 2) != dims[2]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of FFT windows");
			goto fail;
		}
	} else {
		dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
		if(dataF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			goto fail;
		}
	}
	
	#define LAUNCH_PULSAR_ENGINE(IterType) \
        pulsar_engine<IterType>(nStand, nSamps, nFFT, nChan, \
                                (IterType*) PyArray_DATA(data), \
																NULL, \
                                (Complex32*) PyArray_DATA(dataF))
    
    switch( PyArray_TYPE(data) ){
        case( NPY_INT8       ): LAUNCH_PULSAR_ENGINE(int8_t); break;
        case( NPY_COMPLEX64  ): LAUNCH_PULSAR_ENGINE(float);  break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
    
#undef LAUNCH_PULSAR_ENGINE
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(dataF);
	
	return signalsF;
	
fail:
    Py_XDECREF(data);
    Py_XDECREF(dataF);
    
    return NULL;
}

char PulsarEngineRaw_doc[] = PyDoc_STR(\
"Perform a series of Fourier transforms on complex-valued data to get sub-\n\
integration data with linear polarization\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.complex64 (stands by samples) array of data to FFT\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
\n\
Outputs:\n\
 * sub-integration: 2-D numpy.complex64 (stands by channels) of FFT'd data\n\
");


PyObject *PulsarEngineRawWindow(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *window, *signalsF=NULL;
	PyArrayObject *data=NULL, *win=NULL, *dataF=NULL;
	int nChan = 64;
	
	long nStand, nSamps, nFFT;
	
	char const* kwlist[] = {"signals", "window", "LFFT", "signalsF", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO|iO", const_cast<char **>(kwlist), &signals, &window, &nChan, &signalsF)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 2, 2);
	win  = (PyArrayObject *) PyArray_ContiguousFromObject(window,  NPY_DOUBLE,    1, 1);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 2-D complex64");
		goto fail;
	}
	if( win == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input window array to 1-D double");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nSamps = (long) PyArray_DIM(data, 1);
	
	if( PyArray_DIM(win, 0) != nChan ) {
		PyErr_Format(PyExc_RuntimeError, "Window length does not match requested FFT length");
		goto fail;
	}
	
	// Find out how large the output array needs to be and initialize it
	nFFT = nSamps / nChan;
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dims[2] = (npy_intp) (nSamps/nChan);
	if( signalsF != NULL && signalsF != Py_None ) {
		dataF = (PyArrayObject *) PyArray_ContiguousFromObject(signalsF, NPY_COMPLEX64, 3, 3);
		if(dataF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output signalsF array to 3-D complex64");
			Py_XDECREF(data);
			return NULL;
		}
		if(PyArray_DIM(dataF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of stands");
			Py_XDECREF(data);
			Py_XDECREF(dataF);
			return NULL;
		}
		if(PyArray_DIM(dataF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of channels");
			Py_XDECREF(data);
			Py_XDECREF(dataF);
			return NULL;
		}
		if(PyArray_DIM(dataF, 2) != dims[2]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of FFT windows");
			Py_XDECREF(data);
			Py_XDECREF(dataF);
			return NULL;
		}
	} else {
		dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
		if(dataF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			Py_XDECREF(data);
			return NULL;
		}
	}
	
	#define LAUNCH_PULSAR_ENGINE(IterType) \
        pulsar_engine<IterType>(nStand, nSamps, nFFT, nChan, \
                                (IterType*) PyArray_DATA(data), \
																(double*) PyArray_DATA(win), \
                                (Complex32*) PyArray_DATA(dataF))
    
    switch( PyArray_TYPE(data) ){
        case( NPY_INT8       ): LAUNCH_PULSAR_ENGINE(int8_t); break;
        case( NPY_COMPLEX64  ): LAUNCH_PULSAR_ENGINE(float);  break;
        default: PyErr_Format(PyExc_RuntimeError, "Unsupport input data type"); goto fail;
    }
    
#undef LAUNCH_PULSAR_ENGINE
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(win);
	Py_XDECREF(dataF);
	
	return signalsF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(win);
	Py_XDECREF(dataF);
	
	return NULL;
}

char PulsarEngineRawWindow_doc[] = PyDoc_STR(\
"Similar to PulsarEngineRaw but has an additional arugment of the window\n\
function to apply to the data prior to the FFT.\n\
\n\
Input arguments are:\n\
 * signals: 2-D numpy.complex64 (stands by samples) array of data to FFT\n\
 * window: 1-D numpy.float64 array of the window to use (must be the same\n\
   length at the FFT)\n\
\n\
Input keywords are:\n\
 * LFFT: number of FFT channels to make (default=64)\n\
\n\
Outputs:\n\
 * sub-integration: 2-D numpy.complex64 (stands by channels) of FFT'd data\n\
");


PyObject *PhaseRotator(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *signals, *f1, *f2, *signalsF=NULL;
	PyArrayObject *data=NULL, *freq1=NULL, *freq2=NULL, *dataF=NULL;
	double delay;	

	long ij, i, j, k, nStand, nSamps, nChan, nFFT;
	
	char const* kwlist[] = {"signals", "freq1", "freq2", "delays", "signalsF", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOOd|O", const_cast<char **>(kwlist), &signals, &f1, &f2, &delay, &signalsF)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(signals, NPY_COMPLEX64, 3, 3);
	freq1 = (PyArrayObject *) PyArray_ContiguousFromObject(f1, NPY_DOUBLE, 1, 1);
	freq2 = (PyArrayObject *) PyArray_ContiguousFromObject(f2, NPY_DOUBLE, 1, 1);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input signals array to 3-D complex64");
		goto fail;
	}
	if( freq1 == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input f1 array to 1-D double");
		goto fail;
	}
	if( freq2 == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input f2 array to 1-D double");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nChan  = (long) PyArray_DIM(data, 1);
	nFFT   = (long) PyArray_DIM(data, 2);
	nSamps = nChan*nFFT;
	
	// Validate
	if( PyArray_DIM(freq1, 0) != nChan ) {
		PyErr_Format(PyExc_ValueError, "Frequency array 1 has different dimensions than rawSpectra");
		goto fail;
	}
	if( PyArray_DIM(freq2, 0) != nChan ) {
		PyErr_Format(PyExc_ValueError, "Frequency array 2 has different dimensions than rawSpectra");
		goto fail;
	}	
	
	// Find out how large the output array needs to be and initialize it
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dims[2] = (npy_intp) nFFT;
	if( signalsF != NULL && signalsF != Py_None ) {
		dataF = (PyArrayObject *) PyArray_ContiguousFromObject(signalsF, NPY_COMPLEX64, 3, 3);
		if(dataF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output signalsF array to 3-D complex64");
			Py_XDECREF(data);
			return NULL;
		}
		if(PyArray_DIM(dataF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of stands");
			Py_XDECREF(data);
			Py_XDECREF(dataF);
			return NULL;
		}
		if(PyArray_DIM(dataF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of channels");
			Py_XDECREF(data);
			Py_XDECREF(dataF);
			return NULL;
		}
		if(PyArray_DIM(dataF, 2) != dims[2]) {
			PyErr_Format(PyExc_RuntimeError, "signalsF has an unexpected number of FFT windows");
			Py_XDECREF(data);
			Py_XDECREF(dataF);
			return NULL;
		}
	} else {
		dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
		if(dataF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output array");
			goto fail;
		}
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Go!
	long secStart;
	double tempF;
	Complex32 *a, *d;
	double *b, *c;
	a = (Complex32 *) PyArray_DATA(data);
	b = (double *) PyArray_DATA(freq1);
	c = (double *) PyArray_DATA(freq2);
	d = (Complex32 *) PyArray_DATA(dataF);
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(secStart, i, j, k, tempF)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(ij=0; ij<nStand*nChan; ij++) {
			i = ij / nChan;
			j = ij % nChan;
			
			secStart = nSamps*i + nFFT*j;
			if( i/2 == 0 ) {
				tempF = *(b + j);
			} else {
				tempF = *(c + j);
			}
			
			for(k=0; k<nFFT; k++) {
				*(d + secStart + k) = Complex64(*(a + secStart + k)) * exp(TPI*tempF*delay);
			}
		}
	}
	
	Py_END_ALLOW_THREADS
	
	signalsF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(freq1);
	Py_XDECREF(freq2);
	Py_XDECREF(dataF);
	
	return signalsF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(freq1);
	Py_XDECREF(freq2);
	Py_XDECREF(dataF);
	
	return NULL;
}

char PhaseRotator_doc[] = PyDoc_STR(\
"Given the output of PulsarEngineRaw, apply a sub-sample delay as a phase\n\
rotation to each channel\n\
\n\
Input arguments are:\n\
 * signals: 3-D numpy.complex64 (stands by channels by integrations) array\n\
   of data to FFT\n\
 * freq1: 1-D numpy.float64 array of frequencies for each channel for the\n\
   first two stands in signals\n\
 * freq2: 1-D numpy.float64 array of frequencies for each channel for the\n\
   second two stands in signals\n\
 * delay: delay in seconds to apply\n\
\n\
Outputs:\n\
 * signals: 3-D numpy.complex64 (stands by channels by integrations) of the\n\
   phase-rotated spectra data\n\
");
