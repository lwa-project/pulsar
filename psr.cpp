#include "Python.h"
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>
#include <pthread.h>

#ifdef _OPENMP
	#include <omp.h>
	
	// OpenMP scheduling method
	#ifndef OMP_SCHEDULER
	#define OMP_SCHEDULER guided, 16
	#endif
#endif

#define PY_ARRAY_UNIQUE_SYMBOL psr_ARRAY_API
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"

#include "psr.h"


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef SpecMethods[] = {
	{"BindToCore",             (PyCFunction) BindToCore,             METH_VARARGS,               BindToCore_doc             },
	{"BindOpenMPToCores",      (PyCFunction) BindOpenMPToCores,      METH_VARARGS,               BindOpenMPToCores_doc      }, 
	{"PulsarEngineRaw",        (PyCFunction) PulsarEngineRaw,        METH_VARARGS|METH_KEYWORDS, PulsarEngineRaw_doc        },
	{"PulsarEngineRawWindow",  (PyCFunction) PulsarEngineRawWindow,  METH_VARARGS|METH_KEYWORDS, PulsarEngineRawWindow_doc  },
	{"PhaseRotator",           (PyCFunction) PhaseRotator,           METH_VARARGS|METH_KEYWORDS, PhaseRotator_doc           },
	{"ComputeSKMask",          (PyCFunction) ComputeSKMask,          METH_VARARGS,               ComputeSKMask_doc          },
	{"ComputePseudoSKMask",    (PyCFunction) ComputePseudoSKMask,    METH_VARARGS,               ComputePseudoSKMask_doc    },
	{"MultiChannelCD",         (PyCFunction) MultiChannelCD,         METH_VARARGS|METH_KEYWORDS, MultiChannelCD_doc         },
	{"CombineToIntensity",     (PyCFunction) CombineToIntensity,     METH_VARARGS|METH_KEYWORDS, CombineToIntensity_doc     }, 
	{"CombineToLinear",        (PyCFunction) CombineToLinear,        METH_VARARGS|METH_KEYWORDS, CombineToLinear_doc        }, 
	{"CombineToCircular",      (PyCFunction) CombineToCircular,      METH_VARARGS|METH_KEYWORDS, CombineToCircular_doc      }, 
	{"CombineToStokes",        (PyCFunction) CombineToStokes,        METH_VARARGS|METH_KEYWORDS, CombineToStokes_doc        },
	{"OptimizeDataLevels8Bit", (PyCFunction) OptimizeDataLevels8Bit, METH_VARARGS|METH_KEYWORDS, OptimizeDataLevels8Bit_doc },
	{"OptimizeDataLevels4Bit", (PyCFunction) OptimizeDataLevels4Bit, METH_VARARGS|METH_KEYWORDS, OptimizeDataLevels4Bit_doc },
	{NULL,                     NULL,                                 0,                          NULL                       }
};

PyDoc_STRVAR(spec_doc, \
"Extension to take timeseries data and convert it to the frequency domain.\n\
\n\
The functions defined in this module are:\n\
  * PulsarEngineRaw - FFT function for computing a series of Fourier \n\
    transforms for a complex-valued (TBN and DRX) signal from a collection\n\
    of stands/beams all at once.\n\
  * PulsarEngineRawWindow - Similar to PulsarEngineRaw but also requires\n\
    a numpy.float64 array for a window to apply to the data\n\
  * PhaseRotator - Given the output of PulsarEngineRaw, apply a sub-sample\n\
    delay as a phase rotation.\n\
  * ComputeSKMask - Given the output of PulsarEngineRaw compute a mask for\n\
    using spectral kurtosis\n\
  * ComputePseudoSKMask - Similar to ComputeSKMask but for DR spectrometer data\n\
  * MultiChannelCD - Given the output of PulsarEngineRaw apply coherent \n\
    dedispersion to the data\n\
  * CombineToIntensity - Given the output of PulsarEngineRaw compute the total\n\
    intensity for both tunings\n\
  * CombineToLinear - Given the output of PulsarEngineRaw compute XX and YY\n\
    for both tunings\n\
  * CombineToCircular- Given the output of PulsarEngineRaw compute LL and RR\n\
    for both tunings\n\
  * CombineToStokes- Given the output of PulsarEngineRaw compute I, Q, U,\n\
    and V for both tunings\n\
  * OptimizeDataLevels - Given the output of the CombineTo* functions, find\n\
    optimal BZERO and BSCALE values for representing the data as unsigned bytes\n\
    (numpy.uint8)\n\
\n\
See the inidividual functions for more details.");


/*
  Module Setup - Initialization
*/

PyMODINIT_FUNC PyInit__psr(void) {
	char filename[256];
	PyObject *m, *all, *pModule, *pDataPath;

	// Module definitions and functions
	static struct PyModuleDef moduledef = {
		 PyModuleDef_HEAD_INIT, "_psr", spec_doc, -1, SpecMethods
	};
	m  = PyModule_Create(&moduledef);
	if( m == NULL ) {
        return NULL;
    }
	import_array();
	
	// Version information
	PyModule_AddObject(m, "__version__", PyUnicode_FromString("0.6"));
	
	// Function listings
	all = PyList_New(0);
	PyList_Append(all, PyUnicode_FromString("BindToCore"));
	PyList_Append(all, PyUnicode_FromString("BindOpenMPToCores"));
	PyList_Append(all, PyUnicode_FromString("PulsarEngineRaw"));
	PyList_Append(all, PyUnicode_FromString("PulsarEngineRawWindow"));
	PyList_Append(all, PyUnicode_FromString("PhaseRotator"));
	PyList_Append(all, PyUnicode_FromString("ComputeSKMask"));
	PyList_Append(all, PyUnicode_FromString("ComputePseudoSKMask"));
	PyList_Append(all, PyUnicode_FromString("MultiChannelCD"));
	PyList_Append(all, PyUnicode_FromString("CombineToIntensity"));
	PyList_Append(all, PyUnicode_FromString("CombineToLinear"));
	PyList_Append(all, PyUnicode_FromString("CombineToCircular"));
	PyList_Append(all, PyUnicode_FromString("CombineToStokes"));
	PyList_Append(all, PyUnicode_FromString("OptimizeDataLevels8Bit"));
	PyList_Append(all, PyUnicode_FromString("OptimizeDataLevels4Bit"));
	PyList_Append(all, PyUnicode_FromString("useWisdom"));
	PyModule_AddObject(m, "__all__", all);
	
	// LSL FFTW Wisdom
	pModule = PyImport_ImportModule("lsl.common.paths");
	if( pModule != NULL ) {
		pDataPath = PyObject_GetAttrString(pModule, "DATA");
		sprintf(filename, "%s/fftw_wisdom.txt", PyUnicode_AsString(pDataPath));
		read_wisdom(filename, m);
	} else {
		PyErr_Warn(PyExc_RuntimeWarning, "Cannot load the LSL FFTWF wisdom");
	}
	
	return m;
}
