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
  getCoherentSampleSize - Estimate the number of samples needed to 
  successfully apply coherent dedispersion to a data stream.
*/

long getCoherentSampleSize(double centralFreq, double sampleRate, double DM) {
	double delayBand;
	long samples;
	
	delayBand = DM*DCONST * (pow(1e6/(centralFreq-sampleRate/2.0), 2) - pow(1e6/(centralFreq+sampleRate/2.0), 2));
	delayBand *= sampleRate;
	samples = (long) ceil( log(delayBand)/log(2.0) );
	if( samples < 0 ) {
		samples = 0;
	}
	samples = (long) 1<<samples;
	samples *= 2;
	
	return samples;
}


/*
 getFFTChannelFreq - Compute the frequency of a given FFT channel.
*/

double getFFTChannelFreq(long i, long LFFT, double centralFreq, double sampleRate) {
	long N;
	double val;
	
	N = (LFFT-1) / 2 + 1;
	val = sampleRate / (double) LFFT;
	if( i < N ) {
		return i * val + centralFreq;
	} else { 
		return (i - N - LFFT/2) * val + centralFreq;
	}
}


/*
 chirpFunction - Chirp function for coherent dedispersion for a given set of 
 frequencies (in Hz).  Based on Equation (6) of "Pulsar Observations II -- 
 Coherent Dedispersion, Polarimetry, and Timing" By Stairs, I. H.
*/

void chirpFunction(long LFFT, double centralFreq, double sampleRate, double DM, Complex32 *chirp) {
	int i;
	double freqMHz;
	double fMHz0, fMHz1;
	
	// Find the center of the band in MHz
	fMHz0 = centralFreq / 1e6;
	
	// Compute the chirp
	for(i=0; i<LFFT; i++) {
		freqMHz = getFFTChannelFreq(i, LFFT, centralFreq, sampleRate) / 1e6;
		fMHz1 = freqMHz - fMHz0;
		*(chirp + i) = exp(-TPI*DCONST*1e6 * DM*fMHz1*fMHz1 / (fMHz0*fMHz0 * freqMHz));
	}
}


PyObject *MultiChannelCD(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *drxData, *spectraFreq1, *spectraFreq2, *prevData, *nextData, *drxDataF=NULL;
	PyArrayObject *data=NULL, *freq1=NULL, *freq2=NULL;
	PyArrayObject *pData=NULL, *nData=NULL, *dataF=NULL;
	double sRate, DM;
	
	long i, j, k, l, nStand, nChan, nFFT;
	
	char const* kwlist[] = {"rawSpectra", "freq1", "freq2", "sampleRate", "DM", "prevRawSpectra", "nextRawSpectra", "outRawSpectra", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddOO|O", const_cast<char **>(kwlist), &drxData, &spectraFreq1, &spectraFreq2, &sRate, &DM, &prevData, &nextData, &drxDataF)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		goto fail;
	}
	
	// Bring the data into C and make it usable
	data = (PyArrayObject *) PyArray_ContiguousFromObject(drxData, NPY_COMPLEX64, 3, 3);
	freq1 = (PyArrayObject *) PyArray_ContiguousFromObject(spectraFreq1, NPY_DOUBLE, 1, 1);
	freq2 = (PyArrayObject *) PyArray_ContiguousFromObject(spectraFreq2, NPY_DOUBLE, 1, 1);
	pData = (PyArrayObject *) PyArray_ContiguousFromObject(prevData, NPY_COMPLEX64, 3, 3);
	nData = (PyArrayObject *) PyArray_ContiguousFromObject(nextData, NPY_COMPLEX64, 3, 3);
	if( data == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input drxData array to 3-D complex64");
		goto fail;
	}
	if( freq1 == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectraFreq1 to 1-D double");
		goto fail;
	}
	if( freq2 == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input spectraFreq2 to 1-D double");
		goto fail;
	}
	if( pData == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input prevData array to 3-D complex64");
		goto fail;
	}
	if( nData == NULL ) {
		PyErr_Format(PyExc_RuntimeError, "Cannot cast input nextData array to 3-D complex64");
		goto fail;
	}
	
	// Get the properties of the data
	nStand = (long) PyArray_DIM(data, 0);
	nChan  = (long) PyArray_DIM(data, 1);
	nFFT   = (long) PyArray_DIM(data, 2);
	
	// Validate
	if( PyArray_DIM(freq1, 0) != nChan ) {
		PyErr_Format(PyExc_ValueError, "freq1 array has different dimensions than rawSpectra");
		goto fail;
	}
	if( PyArray_DIM(freq2, 0) != nChan ) {
		PyErr_Format(PyExc_ValueError, "freq2 array has different dimensions than rawSpectra");
		goto fail;
	}
	
	if( PyArray_DIM(data, 0) != PyArray_DIM(pData, 0) ) {
		PyErr_Format(PyExc_ValueError, "prevRawSpectra array has different stand dimension than rawSpectra");
		goto fail;
	}
	if( PyArray_DIM(data, 1) != PyArray_DIM(pData, 1) ) {
		PyErr_Format(PyExc_ValueError, "prevRawSpectra array has different channel dimension than rawSpectra");
		goto fail;
	}
	if( PyArray_DIM(data, 2) != PyArray_DIM(pData, 2) ) {
		PyErr_Format(PyExc_ValueError, "prevRawSpectra array has different FFT count dimension than rawSpectra");
		goto fail;
	}
	if( PyArray_DIM(data, 0) != PyArray_DIM(nData, 0) ) {
		PyErr_Format(PyExc_ValueError, "nextRawSpectra array has different stand dimension than rawSpectra");
		goto fail;
	}
	if( PyArray_DIM(data, 1) != PyArray_DIM(nData, 1) ) {
		PyErr_Format(PyExc_ValueError, "nextRawSpectra array has different channel dimension than rawSpectra");
		goto fail;
	}
	if( PyArray_DIM(data, 2) != PyArray_DIM(nData, 2) ) {
		PyErr_Format(PyExc_ValueError, "nextRawSpectra array has different FFT count dimension than rawSpectra");
		goto fail;
	}
	
	// Find out how large the output arrays needs to be and initialize then
	npy_intp dims[3];
	dims[0] = (npy_intp) nStand;
	dims[1] = (npy_intp) nChan;
	dims[2] = (npy_intp) nFFT;
	if( drxDataF != NULL && drxDataF != Py_None ) {
		dataF = (PyArrayObject *) PyArray_ContiguousFromObject(drxDataF, NPY_COMPLEX64, 3, 3);
		if(dataF == NULL) {
			PyErr_Format(PyExc_RuntimeError, "Cannot cast output outRawSpectra array to 3-D complex64");
			goto fail;
		}
		if(PyArray_DIM(dataF, 0) != dims[0]) {
			PyErr_Format(PyExc_RuntimeError, "outRawSpectra has an unexpected number of stands");
			goto fail;
		}
		if(PyArray_DIM(dataF, 1) != dims[1]) {
			PyErr_Format(PyExc_RuntimeError, "outRawSpectra has an unexpected number of channels");
			goto fail;
		}
		if(PyArray_DIM(dataF, 2) != dims[2]) {
			PyErr_Format(PyExc_RuntimeError, "outRawSpectra has an unexpected number of FFT windows");
			goto fail;
		}
	} else {
		dataF = (PyArrayObject*) PyArray_ZEROS(3, dims, NPY_COMPLEX64, 0);
		if(dataF == NULL) {
			PyErr_Format(PyExc_MemoryError, "Cannot create output data array");
			goto fail;
		}
	}
	
	Py_BEGIN_ALLOW_THREADS
	
	// Get access to the frequency information
	double *cf1, *cf2;
	cf1 = (double *) PyArray_DATA(freq1);
	cf2 = (double *) PyArray_DATA(freq2);
	
	// Create the FFTW plans
	long N;
	Complex32 *inP;
	fftwf_plan *plansF, *plansB;
	plansF = (fftwf_plan *) malloc(nStand/2*nChan*sizeof(fftwf_plan));
	plansB = (fftwf_plan *) malloc(nStand/2*nChan*sizeof(fftwf_plan));
	for(j=0; j<nChan; j++) {
		for(i=0; i<nStand; i+=2) {
			//Compute the number of FFT channels to use
			if( i/2 == 0 ) {
				N = getCoherentSampleSize(*(cf1 + j), sRate, DM);
			} else {
				N = getCoherentSampleSize(*(cf2 + j), sRate, DM);
			}
			
			inP = (Complex32 *) fftwf_malloc(N*sizeof(Complex32));
			*(plansF + i/2*nChan + j) = fftwf_plan_dft_1d(N,
																										reinterpret_cast<fftwf_complex*>(inP),
																										reinterpret_cast<fftwf_complex*>(inP),
																										FFTW_FORWARD, FFTW_ESTIMATE);
			*(plansB + i/2*nChan + j) = fftwf_plan_dft_1d(N,
																										reinterpret_cast<fftwf_complex*>(inP),
																										reinterpret_cast<fftwf_complex*>(inP),
																										FFTW_BACKWARD, FFTW_ESTIMATE);
			fftwf_free(inP);
		}
	}
	
	// Go!
	long nSets, start, stop, secStartX, secStartY;
	double cFreq;
	Complex32 *chirp;
	Complex32 *d0, *d1, *d2, *dF;
	
	Complex32 *inX, *inY;
	
	d0 = (Complex32 *) PyArray_DATA(pData);
	d1 = (Complex32 *) PyArray_DATA(data);
	d2 = (Complex32 *) PyArray_DATA(nData);
	dF = (Complex32 *) PyArray_DATA(dataF);
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(secStartX, secStartY, i, j, k, l, N, nSets, start, stop, cFreq, chirp, inX, inY)
	#endif
	{
		#ifdef _OPENMP
			#pragma omp for schedule(OMP_SCHEDULER)
		#endif
		for(j=0; j<nChan; j++) {
			for(i=0; i<nStand; i+=2) {
				// Section start offset
				secStartX = i*nChan*nFFT     + j*nFFT;
				secStartY = (i+1)*nChan*nFFT + j*nFFT;
				
				// Get the correct center frequency to use
				if( i/2 == 0 ) {
					cFreq = *(cf1 + j);
				} else {
					cFreq = *(cf2 + j);
				}
				
				//Compute the number of FFT channels to use
				N = getCoherentSampleSize(cFreq, sRate, DM);
				
				// Compute the number of windows we need to use for CD
				nSets = nFFT / N;
				
				// Compute the chirp function
				chirp = (Complex32 *) malloc(N*sizeof(Complex32));
				chirpFunction(N, cFreq, sRate, DM, chirp);
				
				// Create the FFTW array
				inX = (Complex32 *) fftwf_malloc(N*sizeof(Complex32));
				inY = (Complex32 *) fftwf_malloc(N*sizeof(Complex32));
				
				// Loop over the sets
				for(l=0; l<2*nSets+1; l++) {
					start = N/2*l - N/4;
					stop = start + N;
					
					// Load in the data
					if( start < 0 ) {
						// "Previous" buffering
						for(k=0; k<-start; k++) {
							inX[k] = *(d0 + secStartX + k + nFFT + start);
							inY[k] = *(d0 + secStartY + k + nFFT + start);
						}
						
						// Current data
						for(k=-start; k<N; k++) {
							inX[k] = *(d1 + secStartX + k + start);
							inY[k] = *(d1 + secStartY + k + start);
						}
						
					} else if( stop > nFFT ) {
						// Current data
						if( start < nFFT ) {
							for(k=0; k<nFFT-start; k++) {
								inX[k] = *(d1 + secStartX + k + start);
								inY[k] = *(d1 + secStartY + k + start);
							}
						}
						
						// "Next" buffering
						for(k=nFFT-start; k<N; k++) {
							inX[k] = *(d2 + secStartX + k - nFFT + start);
							inY[k] = *(d2 + secStartY + k - nFFT + start);
						}
						
					} else {
						// Current data
						for(k=0; k<N; k++) {
							inX[k] = *(d1 + secStartX + k + start);
							inY[k] = *(d1 + secStartY + k + start);
						}
						
					}
					
					// Forward FFT
					fftwf_execute_dft(*(plansF + i/2*nChan + j),
														reinterpret_cast<fftwf_complex*>(inX),
														reinterpret_cast<fftwf_complex*>(inX));
					fftwf_execute_dft(*(plansF + i/2*nChan + j),
													  reinterpret_cast<fftwf_complex*>(inY),
														reinterpret_cast<fftwf_complex*>(inY));
					
					// Chirp
					for(k=0; k<N; k++) {
						inX[k] *= *(chirp + k) / (float) N;
						inY[k] *= *(chirp + k) / (float) N;
					}
					
					// Backward FFT
					fftwf_execute_dft(*(plansB + i/2*nChan + j),
					                  reinterpret_cast<fftwf_complex*>(inX),
														reinterpret_cast<fftwf_complex*>(inX));
					fftwf_execute_dft(*(plansB + i/2*nChan + j),
													  reinterpret_cast<fftwf_complex*>(inY),
														reinterpret_cast<fftwf_complex*>(inY));
					
					// Save
					start = N/2*l;
					stop = start + N/2;
					if( stop > nFFT ) {
						stop = nFFT;
					}
					
					for(k=start; k<stop; k++) {
						*(dF + secStartX + k) = inX[k - start + N/4];
						*(dF + secStartY + k) = inY[k - start + N/4];
					}
					
					if( stop == nFFT ) {
						break;
					}
					
				}
				
				// Cleanup
				free(chirp);
				fftwf_free(inX);
				fftwf_free(inY);
			}
		}
	}
	
	// Cleanup
	for(j=0; j<nChan; j++) {
		for(i=0; i<nStand; i+=2) {
			fftwf_destroy_plan(*(plansF + i/2*nChan + j));
			fftwf_destroy_plan(*(plansB + i/2*nChan + j));
		}
	}
	free(plansF);
	free(plansB);
	
	Py_END_ALLOW_THREADS
	
	drxDataF = Py_BuildValue("O", PyArray_Return(dataF));
	
	Py_XDECREF(data);
	Py_XDECREF(freq1);
	Py_XDECREF(freq2);
	Py_XDECREF(pData);
	Py_XDECREF(nData);
	Py_XDECREF(dataF);
	
	return drxDataF;
	
fail:
	Py_XDECREF(data);
	Py_XDECREF(freq1);
	Py_XDECREF(freq2);
	Py_XDECREF(pData);
	Py_XDECREF(nData);
	Py_XDECREF(dataF);
	
	return NULL;
}

char MultiChannelCD_doc[] = PyDoc_STR(\
"Given the output of one of the 'PulsarEngine' functions and information about\n\
the time and frequencies, apply coherent dedispersion to the data and return a\n\
two-element tuple giving the time and dedispersed data.\n\
\n\
Input arguments are:\n\
  * rawSpectra - 3-D numpy.complex64 (stands by channels by samples) of raw\n\
    spectra data generated by 'PulsarEngineRaw' or 'PulsarEngineWindow'\n\
  * freq1 - 1-D numpy.float64 (channels) array of frequencies for the first\n\
    set of two stands in Hz\n\
  * freq2 - 1-D numpy.float64 (channels) array of frequency for the second\n\
    set of two stands in Hz\n\
  * sampleRate - Channelized data sample rate in Hz\n\
  * DM - dispersion measure in pc cm^-3 to dedisperse at\n\
  * prevRawSpectra - 3-D numpy.complex64 (stands by channels by samples) of\n\
    the previously input data\n\
  * prevRawSpectra - 3-D numpy.complex64 (stands by channels by samples) of\n\
    the following input data\n\
\n\
Outputs:\n\
  * dedispRawSpectra - 3-D numpy.complex64 (stands by channels by samples) of\n\
    the coherently dedispersed spectra\n\
\n\
.. note::\n\
\tThere are a few things that look a little strange here.  First, the\n\
\tsample rate specified is that of the output spectra.  This is related\n\
\tto the input data range and channel count via <input rate> / <channel\n\
\tcount>.  Second, the time/prevTime/nextTime, etc. variables might be\n\
\ta little confusing.  Consider the following time/data flow in a file:\n\
\t  t0/d0, t1/d1, t2/d2, t3/d3, ...\n\
\tFor the first call to MultiDisp() the arguments are:\n\
\t  * t1/d1 are time/rawSpectra\n\
\t  * t0/d0 are prevTime/prevRawSpectra\n\
\t  * t2/d2 are nextTime/nextRawSpectra\n\
\tand for the subsequent call to MultiDisp() they are:\n\
\t  * t2/d2 are time/rawSpectra\n\
\t  * t1/d1 are prevTime/prevRawSpectra\n\
\t  * t3/d3 are nextTime/nextRawSpectra\n\
");
