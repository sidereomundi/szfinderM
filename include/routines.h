//$Id: routines.h,v 1.4 2012/11/27 13:50:15 jiayi Exp $
// Declare the routines in used
#ifndef _ROUTINES_H_
#define _ROUTINES_H_

#include "fftw3.h"
struct simmap;

// math related
#define Squ(x) ((x)*(x))
#define Cube(x) ((x)*(x)*(x))
#define Quad(x) ((x)*(x)*(x)*(x))

// sub_setup.c
void initializemap(struct simmap*);   // allocate the simulation map
void cleanmap(struct simmap*);       // delocate the simulation map
float *createmap(struct simmap*);    // create a map with same dimensions as the simulation map
void savemap(struct simmap*, char*); // save map image to fits
void mapping(struct simmap*, float*);

// sub_clinput.c
void clinput(float**, float **, float**, float*, float*, int*);
        // !!allocate and input the cl data!!

// sub_mkdirtmap.c
void mkdirtmap(struct simmap*, float*, fftw_complex*);  // add noise and beam effect
void removemean(struct simmap*);

// sub_filter.c
void filter(fftw_complex*, fftw_complex*, int, int, float*);         // filter two maps (origin map, filter map) return to another map

// sub_fft.c
void fftread(long, fftw_complex*, char*);       // read fft written by fftsave()
void fft(float*, int, int, fftw_complex*);      // convert real image to complex 2D array (half complex plane)
void fftc(float*, int, int, fftw_complex*);     // convert real image to complex 2D array (full complex plane)

// sub_fftsave.c
void fftsave(fftw_complex*, long, char *);
// sub_fft2image.c
void fft2image(fftw_complex*, int, int, float*);  // convert to real image based on half complex plane
void fftc2image(fftw_complex*, int, int, char*);  // convert to real image based on full complex plane and output to image
void fft2imageaddnoise(fftw_complex *,int,int,float *,float *); // convert to real image based on

// sub_stdroutines.c
void saveimagef(char*, float*, int, int);     // write fits file
void readimage(char*, float**, int*, int*);     // read fits file 
        // !! allocate the image array !!
double thetasize(int);    // return the theta_core radius for spt filter [arcmin]

// sub_mkcmb.c
void mkcmb(int, long*, struct simmap*, float*, float*, float*, float*, float, float); // make CMB temperature map

// sub_mkbolo.c
void mkbolo(struct simmap*, double, float*, float*, double, long*);   // convert sz-y, noise to their bolometric unit (temperature as Tcmb)

// sub_tauvalue.c
double tauvalue(int, int, double);  // return cluster's sz profile

#endif
