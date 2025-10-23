// $Id: main_mkoneszmap.c,v 1.2 2013/05/10 15:14:31 jiayi Exp $
// make mock sz map

#include "fftw3.h"
#include <stdio.h>
#include <stdlib.h>

#include "routines.h"
#include "parameters.h"
#include "structures.h"

/**
 * \brief make simulated flat sky
 * 
 * usage:
 *        mkoneszmap <rseed> <frequency> <white_noise> <beam> <Ymap> <output> <rseed noise> 
 * input:
 *        rseed:
 *        frequency:
 *        white_noise:
 *        beam:
 *        Ymap:
 *        rseed_noise:
 * default input:
 *        
*/
int main(int argc, char *argv[])
{
  long rseed,rseednoise;
  struct simmap imap;
  float *szeimage;
  int xdim,ydim,ntot;
  char filename[500],beamname[500];
  double freq,wnoise;

  float *cl_l,*cl_p,*cl_py2,yp1,ypn;
  int cllength;
  float *cmbimage;
  fftw_complex *fbeam;
  float *noise;
   
// initial set up
  if(argc != 8){
    printf("mkoneszmap <rseed> <frequency> <white_noise> <beam file> <Y map file> <output> <rseed noise>\n");
    exit(1);
  }

  sscanf(argv[1],"%ld",&rseed);
  sscanf(argv[2],"%lf",&freq);
  sscanf(argv[3],"%lf",&wnoise);
  sscanf(argv[4],"%s",beamname);  // read in beamname_fft.dat
  sscanf(argv[5],"%s",filename);
  sscanf(argv[7],"%ld",&rseednoise);
  rseed = -rseed;
  rseednoise = -rseednoise;
  
  initializemap(&imap);
  // SZE
  printf("# Read the SZE image %s\n",filename);
  readimage(filename,&szeimage,&xdim,&ydim);
  if((xdim!=NDIM)||(ydim!=NDIM)){
    printf("# size of image is wrong!\n");
    exit(1);
  }
  // put the sze image into imap (could add boundary layer)
  mapping(&imap,szeimage);
  xdim = imap.nsidepix;
  ydim = imap.nsidepix;

  // CMB
  printf("# Input CMB \n");
  clinput(&cl_l,&cl_p,&cl_py2,&yp1,&ypn,&cllength);
  cmbimage = (float *)malloc(sizeof(float)*imap.npix);
  mkcmb(cllength,&rseed,&imap,cmbimage,cl_l,cl_p,cl_py2,yp1,ypn);
  saveimagef("cmb.fits",cmbimage,xdim,ydim);

  // mkbolo
  printf("# Bolometry @ %lf GHz\n",freq);
  noise=(float*)malloc(sizeof(float)*imap.npix);
  mkbolo(&imap,freq,cmbimage,noise,wnoise,&rseednoise);
  saveimagef("noise.fits",noise,xdim,ydim);
  // beam effect
  printf("# Adding beam and white noise\n");
  ntot = ydim*(xdim/2+1);
  fbeam=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntot);
  sprintf(filename,"%s_fft.dat",beamname);
  fftread(ntot,fbeam,filename);
  mkdirtmap(&imap,noise,fbeam);

  // substract the Tcmb
  removemean(&imap);

  // Output
  printf("# Output to %s\n",argv[6]);
  savemap(&imap,argv[6]);

  // clean
  free(cmbimage);
  free(szeimage);
  free(noise);
  fftw_free(fbeam);
  free(cl_l);
  free(cl_p);
  free(cl_py2);
  return 0;
}
