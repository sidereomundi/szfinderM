//<<<<<<< parameters.h
//$Id: parameters.h,v 1.10 2012/11/27 13:50:15 jiayi Exp $
//=======
//$Id: parameters.h,v 1.13 2013/10/16 10:32:02 jiayi Exp $
//>>>>>>> 1.13
// Define various parameters for the simulator

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

// general parameters
#define clfile "./cl.dat"  // the power spectrum of CMB

// mock related
#define NDIM 4800                 // pixel on one side

#define FIELDSIZE 20              // field size [degree]
#define EXTENDEDDIM 0           // extended pixel number per side

#define FREQUENCY150 154.477        // 150 GHz band frequency [GHz]
#define DT150 18.0                // 150 GHz band noise [K/arcmin]
#define BEAM150FILE "beam150T"  // 150 GHz band beam (Gaussian 1.1 arcmin)
#define FREQUENCY90 97.6473        // 95 GHz band frequency [GHz]
#define DT90 44.0                // 95 GHz band noise [K/arcmin]
#define BEAM90FILE "beam95T"  // 95 GHz band beam (Gaussian 1.6 arcmin)

// cluster finder related
#define FILTERSIGFILE "sigfilter.dat"  // variance of the filter
#define FILTERMSIGFILE "sigfilterM.dat"  // variance of the filter (multifrequency)
#define NFILTER 12                      // number of filter 
#define FILTERFILE "filter"        // file name head for filter  // chagned from filter150_
#define BETA 1                         // parameter for beta profile
#define THETA_START 0.25                // theta_core for first beta profile (arcmin)
#define THETA_STEP  0.25                // theta_core steps for beta profile
#define CLUSTERLIST "clusterlist.dat"  // cluster list formatted from the simulation
#define SIGMA_THRESHOD 5.              // cluster sn limit
#define SZPEAKLIST "szpeaks.dat"       // output sz peak list
// physics parameters
#define T_CMB 2.728		// Temperature of CMB [K]
#define KB      1.38066e-16	// Boltzmann constant in CGS
#define H	6.6260755e-27	// Planck constant in CGS
#define C_CGS   2.9979247e+10 	/* speed of light in cm/s */

// math related
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// temporary definitions

//#define pix2arcmin 0.08496094 // 2.9*60/2048 arcmin
#endif
