#include "Python.h"
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <pthread.h>

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
  Core binding function - based off the corresponding bifrost function
*/

#if defined __APPLE__ && __APPLE__

// Based on information from:
//   http://www.hybridkernel.com/2015/01/18/binding_threads_to_cores_osx.html

#include <sys/types.h>
#include <sys/sysctl.h>
#include <mach/mach_init.h>
#include <mach/thread_policy.h>
#include <mach/thread_act.h>

typedef struct cpu_set {
  uint32_t    count;
} cpu_set_t;

static inline void
CPU_ZERO(cpu_set_t *cs) { cs->count = 0; }

static inline void
CPU_SET(int num, cpu_set_t *cs) { cs->count |= (1 << num); }

static inline int
CPU_ISSET(int num, cpu_set_t *cs) { return (cs->count & (1 << num)); }

static inline int
CPU_COUNT(cpu_set_t *cs) {
	int count = 0;
	for(int i=0; i<8*sizeof(cpu_set_t); i++) {
		count += CPU_ISSET(i, cs) ? 1 : 0;
	}
	return count;
}

int pthread_getaffinity_np(pthread_t thread,
	                         size_t    cpu_size,
                           cpu_set_t *cpu_set) {
  thread_port_t mach_thread;
	mach_msg_type_number_t count = THREAD_AFFINITY_POLICY_COUNT;
	boolean_t get_default = false;
	
	thread_affinity_policy_data_t policy;
	mach_thread = pthread_mach_thread_np(thread);
	thread_policy_get(mach_thread, THREAD_AFFINITY_POLICY,
                    (thread_policy_t)&policy, &count,
									  &get_default);
	cpu_set->count |= (1<<(policy.affinity_tag));
	return 0;
}

int pthread_setaffinity_np(pthread_t thread,
	                         size_t    cpu_size,
                           cpu_set_t *cpu_set) {
  thread_port_t mach_thread;
  int core = 0;

  for (core=0; core<8*cpu_size; core++) {
    if (CPU_ISSET(core, cpu_set)) break;
  }
  thread_affinity_policy_data_t policy = { core };
  mach_thread = pthread_mach_thread_np(thread);
  thread_policy_set(mach_thread, THREAD_AFFINITY_POLICY,
                    (thread_policy_t)&policy, 1);
  return 0;
}

#endif

/*
 setCore - Internal function to bind a thread to a particular core, or unbind
 it if core is set to -1.  Returns 0 is successful, non-zero otherwise.
*/         

int setCore(int core) {
#if (defined __linux__ && __linux__) || (defined __APPLE__ && __APPLE__)
	int ncore;
	cpu_set_t cpuset;
	ncore = sysconf(_SC_NPROCESSORS_ONLN);
	
	// Basic validation
	if( core >= ncore || (core < 0 && core != -1) ) {
		return -100;
	}
	
	CPU_ZERO(&cpuset);
	if( core >= 0 ) {
		CPU_SET(core, &cpuset);
	} else {
		for(core=0; core<ncore; ++core) {
			CPU_SET(core, &cpuset);
		}
	}
	
	pthread_t tid = pthread_self();
	return pthread_setaffinity_np(tid, sizeof(cpu_set_t), &cpuset);
#else
	return -101;
#endif
}

/*
 getCore - Internal function to get the core binding of a thread.  Returns 
 0 is successful, non-zero otherwise.
*/

int getCore(int* core) {
#if (defined __linux__ && __linux__) || (defined __APPLE__ && __APPLE__)
	int ret, c, ncore;
	cpu_set_t cpuset;
	ncore = sysconf(_SC_NPROCESSORS_ONLN);
	
	pthread_t tid = pthread_self();
	ret = pthread_getaffinity_np(tid, sizeof(cpu_set_t), &cpuset);
	if( ret == 0 ) {
		if( CPU_COUNT(&cpuset) > 1 ) {
			*core = -1;
			return 0;
		} else {
			for(c=0; c<ncore; ++c ) {
				if(CPU_ISSET(c, &cpuset)) {
					*core = c;
					return 0;
				}
			}
		}
	}
	return ret;
#else
	return -101;
#endif
}

/* getCoreCount - Internal function to return the number of core available.  Returns >=1
   if successful, <1 otherwise.
*/

int getCoreCount(void) {
#if (defined __linux__ && __linux__) || (defined __APPLE__ && __APPLE__)
	return sysconf(_SC_NPROCESSORS_ONLN);
#else
	return -101;
#endif
}

PyObject *BindToCore(PyObject *self, PyObject *args, PyObject *kwds) {
	int ret, core = -1;
	
	if(!PyArg_ParseTuple(args, "i", &core)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	ret = setCore(core);
	
	if(ret == 0) {
		Py_RETURN_TRUE;
	} else if(ret == -100) {
		PyErr_Format(PyExc_ValueError, "Invalid core: %d", core);
		return NULL;
	} else if(ret == -101) {
		PyErr_Warn(PyExc_RuntimeWarning, "Changing of the thread core binding is not supported on this OS");
		Py_RETURN_FALSE;
	} else {
		PyErr_Format(PyExc_RuntimeError, "Cannot change core binding to %d", core);
		return NULL;
	}
}

char BindToCore_doc[] = PyDoc_STR(\
"Bind the current thread to the specified core.\n\
\n\
Input arguments are:\n\
 * core: scalar int core to bind to\n\
\n\
Outputs:\n\
 * True, if successful\n\
");


PyObject *BindOpenMPToCores(PyObject *self, PyObject *args, PyObject *kwds) {
	PyObject *cores, *core;
	int ret, nthread, t, tid, ncore, old_core, c;
	
	if(!PyArg_ParseTuple(args, "O", &cores)) {
		PyErr_Format(PyExc_RuntimeError, "Invalid parameters");
		return NULL;
	}
	
	if(!PyList_Check(cores)) {
		PyErr_Format(PyExc_TypeError, "Invalid parameters");
		return NULL;
	}
	
	nthread = (int) PyList_Size(cores);
	ncore = getCoreCount();
	if( ncore == -101 ) {
		PyErr_Warn(PyExc_RuntimeWarning, "Changing of the thread core binding is not supported on this OS");
		Py_RETURN_FALSE;
	}
	for(t=0; t<nthread; t++) {
		core = PyList_GetItem(cores, t);
		if(!PyInt_Check(core)) {
			PyErr_Format(PyExc_TypeError, "Invalid parameters");
			return NULL;
		}
		c = (int) PyInt_AsLong(core);
		if( c >= ncore || (c < 0 && c != -1 ) ) {
			PyErr_Format(PyExc_ValueError, "Invalid core for thread %d: %d", t+1, c);
			return NULL;
		}
	}
	
	ret = getCore(&old_core);
	if(ret == -101) {
		PyErr_Warn(PyExc_RuntimeWarning, "Changing of the thread core binding is not supported on this OS");
		Py_RETURN_FALSE;
	} else if(ret != 0) {
		PyErr_Format(PyExc_RuntimeError, "Cannot get core binding");
		return NULL;
	}
	
	ret = setCore(-1);
	if(ret != 0) {
		PyErr_Format(PyExc_RuntimeError, "Cannot unbind current thread");
		return NULL;
	}
	
	#ifdef _OPENMP
		#pragma omp parallel default(shared) private(tid, core, c)
		omp_set_num_threads(nthread);
		
		#pragma omp parallel for schedule(static, 1)
		for(t=0; t<nthread; ++t) {
			tid = omp_get_thread_num();
			core = PyList_GetItem(cores, tid);
			c = (int) PyInt_AsLong(core);
			ret |= setCore(c);
		}
	#else
		ret = -101;
	#endif
	if(ret == -101) {
		PyErr_Warn(PyExc_RuntimeWarning, "Changing of the thread core binding is not supported on this OS");
		Py_RETURN_FALSE;
	} else if(ret != 0) {
		PyErr_Format(PyExc_RuntimeError, "Cannot set all OpenMP thread core bindings");
		return NULL;
	}
	
	if(old_core != -1) {
		ret = setCore(old_core);
		if(ret != 0) {
			PyErr_Format(PyExc_RuntimeError, "Cannot re-bind current thread to %d", old_core);
			return NULL;
		}
	}
	
	Py_RETURN_TRUE;
}

char BindOpenMPToCores_doc[] = PyDoc_STR(\
"Bind OpenMP threads to the provided list of cores.\n\
\n\
Input arguments are:\n\
 * cores: list of int cores for OpenMP thread binding\n\
\n\
Outputs:\n\
  * True, if successful\n\
");
