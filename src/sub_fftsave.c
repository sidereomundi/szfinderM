// $Id: sub_fftsave.c,v 1.1 2012/06/04 14:03:20 jiayi Exp $
// save the result from fft() (half complex plane)
#include "fftw3.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void fftsave(fftw_complex *fdata, long ntot, char *fname)
{
  FILE *fft;
  char datafft[200];
  int i;
  /* Set up the file name */
  if (strncmp(&(fname[strlen(fname)-4]),".dat",4))
    sprintf(datafft,"%s_fft.dat",fname);
  else {
    for(i=0;i<strlen(fname)-5;i++)
      datafft[i]=fname[i];
    datafft[i]='\0';
    strcat(datafft,"_fft.dat");
  }
  printf(" * saving fft result\n");
  fft=fopen(datafft,"w");
  
  /* unformatted output	*/
  fwrite(fdata,sizeof(fftw_complex),ntot,fft);
  fclose(fft);
}
