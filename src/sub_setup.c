//$Id: sub_setup.c,v 1.3 2012/06/07 16:24:40 jiayi Exp $
// Setup for various structure

#include <stdio.h>
#include <stdlib.h>

#include "parameters.h"
#include "structures.h"
#include "routines.h"

void initializemap(struct simmap *pmap)
{
  pmap->pixsize = (double)FIELDSIZE/NDIM * 60;  // in arcmin
  printf("# pixel size is %f [arcmin]\n",pmap->pixsize);
  pmap->nsidepix = NDIM+EXTENDEDDIM*2;
  pmap->npix = pmap->nsidepix*pmap->nsidepix;
  pmap->img = (float*)malloc(sizeof(float)*pmap->npix);
  if (pmap->img == NULL) {
    printf("!! fail to allocate memory for simmap->img !!\n");
    exit(1);
  }
}

float *createmap(struct simmap *pmap)
{
  float *temp;
  temp = (float*)malloc(sizeof(float)*pmap->npix);
  if (temp == NULL) {
    printf("!! fail to allocate memory for createmap !!\n");
    exit(1);
  }
  return temp;
}

void savemap(struct simmap *pmap, char *fname)
{
  saveimagef(fname,pmap->img,pmap->nsidepix,pmap->nsidepix);
}

void cleanmap(struct simmap *pmap)
{
  pmap->pixsize = 0;
  pmap->nsidepix = 0;
  pmap->npix = 0;
  free(pmap->img);
  pmap->img = NULL;
}

void mapping(struct simmap *pmap, float *szimg)
{
  long ip,ix,iy;
  long sx,sy,sp;
  for (ix=0;ix<pmap->nsidepix;ix++){
    sx = ix-EXTENDEDDIM;
    if ((sx<0)||(sx>=NDIM))
      for(iy=0;iy<pmap->nsidepix;iy++){
	ip = ix + iy*pmap->nsidepix;
	pmap->img[ip]=0;
      }
    else
      for(iy=0;iy<pmap->nsidepix;iy++){
	ip = ix + iy*pmap->nsidepix;
	sy = iy-EXTENDEDDIM;
	if ((sy<0)||(sy>=NDIM))
	  pmap->img[ip]=0;
	else{
	  sp = sx + sy*NDIM;
	  pmap->img[ip]=szimg[sp];
	}
      }
  }
}
