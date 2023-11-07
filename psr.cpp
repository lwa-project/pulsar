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

#include "psr.hpp"


/*
  Module Setup - Function Definitions and Documentation
*/

static PyMethodDef psr_methods[] = {
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

PyDoc_STRVAR(psr_doc, \
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

static int psr_exec(PyObject *module) {
		import_array();
		
		// Version information
		PyModule_AddObject(module, "__version__", PyUnicode_FromString("0.6"));
		
		// Function listings
		PyObject* all = PyList_New(0);
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
		PyModule_AddObject(module, "__all__", all);
		
		// LSL FFTW Wisdom
		PyObject* pModule = PyImport_ImportModule("lsl.common.paths");
		if( pModule != NULL ) {
				PyObject* pDataPath = PyObject_GetAttrString(pModule, "DATA");
				if( pDataPath != NULL ) {
						char filename[256];
						sprintf(filename, "%s/fftw_wisdom.txt", PyUnicode_AsString(pDataPath));
						read_wisdom(filename, module);
				}
				Py_XDECREF(pDataPath);
		} else {
				PyErr_Warn(PyExc_RuntimeWarning, "Cannot load the LSL FFTWF wisdom");
		}
		Py_XDECREF(pModule);
		return 0;
}

static PyModuleDef_Slot psr_slots[] = {
    {Py_mod_exec, (void *)&psr_exec},
    {0,           NULL}
};

static PyModuleDef psr_def = {
    PyModuleDef_HEAD_INIT,    /* m_base */
    "_psr",                   /* m_name */
    psr_doc,                  /* m_doc */
    0,                        /* m_size */
    psr_methods,              /* m_methods */
    psr_slots,                /* m_slots */
    NULL,                     /* m_traverse */
    NULL,                     /* m_clear */
    NULL,                     /* m_free */
};

PyMODINIT_FUNC PyInit__psr(void) {
	return PyModuleDef_Init(&psr_def);
}
