// $Id: main_mkszmap.c,v 1.6 2013/02/12 14:23:15 jiayi Exp $
// make mock sz map

#include "fftw3.h"
#include <stdio.h>
#include <stdlib.h>

#include "routines.h"
#include "parameters.h"
#include "structures.h"

int main(int argc, char *argv[])
{
  long rseed;
  struct simmap imap;
  float *szeimage;
  int xdim,ydim,ntot;
  char filename[500];
  
  float *cl_l,*cl_p,*cl_py2,yp1,ypn;
  int cllength;
  float *cmbimage;
  fftw_complex *fbeam;
  float *noise;
	
// initial set up
  if(argc<2){
    printf("mkszmap <rseed> <filename of sz-Y map>\n");
    exit(1);
  }
  else if(argc<3){
    printf("using default file: sznull.fits");
    sprintf(filename,"sznull.fits");
  }
  else sscanf(argv[2],"%s",filename);
  sscanf(argv[1],"%ld",&rseed);
  rseed = -rseed;
  
  initializemap(&imap);
  // SZE
  printf("# Read the SZE image %s\n",filename);
  readimage(filename,&szeimage,&xdim,&ydim);
  if((xdim!=NDIM)||(ydim!=NDIM)){
    printf("# size of image is wrong!\n");
    exit(1);
  }
  mapping(&imap,szeimage);
  xdim = imap.nsidepix;
  ydim = imap.nsidepix;
  //  savemap(&imap,"szremapping.fits");  remapping is correct

  // CMB
  printf("# Input CMB \n");
  clinput(&cl_l,&cl_p,&cl_py2,&yp1,&ypn,&cllength);
  cmbimage = (float *)malloc(sizeof(float)*imap.npix);
  mkcmb(cllength,&rseed,&imap,cmbimage,cl_l,cl_p,cl_py2,yp1,ypn);
  //  saveimagef("cmb.fits",cmbimage,xdim,ydim);
  // 150 GHz
  // mkbolo
  printf("# Bolometry\n");
  noise=(float*)malloc(sizeof(float)*imap.npix);
  mkbolo(&imap,FREQUENCY150,cmbimage,noise,DT150,&rseed);
  //  saveimagef("NoNoise_beam.fits",imap.img,xdim,ydim);
  //  saveimagef("noise.fits",noise,xdim,ydim);
  // beam effect
  printf("# Adding beam and white noise\n");
  ntot = ydim*(xdim/2+1);
  fbeam=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*ntot);
  sprintf(filename,"%s_fft.dat",BEAM150FILE);
  fftread(ntot,fbeam,filename);
  mkdirtmap(&imap,noise,fbeam);
  removemean(&imap);      
  // Output
  printf("# Output \n");
  //	sprintf(filename,"pure_noise.fits");  // for sz-null image
  sprintf(filename,"spt150.fits");
  savemap(&imap,filename);

  // for 90 GHz
  printf("Addition for 90 GHz\n");
  mapping(&imap,szeimage);
  mkbolo(&imap,FREQUENCY90,cmbimage,noise,DT90,&rseed);
  sprintf(filename,"%s_fft.dat",BEAM90FILE);
  fftread(ntot,fbeam,filename);
  mkdirtmap(&imap,noise,fbeam);
  removemean(&imap);
  sprintf(filename,"spt95.fits");
  savemap(&imap,filename);

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
