// test file


#include <stdio.h>
#include <stdlib.h>
#include "nr.h"
#include "routines.h"
#include "structures.h"
#include "parameters.h"

int main()
{
  // test on structure
  /*  struct simmap imap;
  initializemap(&imap);
  printf("map size = %ld\n",imap.npix);
  finalizemap(&imap);
  */
  // test on clinput  --DONE--
  /*  float *cl_l,*cl_p,*cl_py2,yp1,ypn;
  int cllength;
  clinput(&cl_l,&cl_p,&cl_py2,&yp1,&ypn,&cllength);
  free(cl_l);
  free(cl_p);
  free(cl_py2);
  */
  // test on fits library  // by-product: test on screen color
  /*  float *img;
  img = (float*)malloc(sizeof(float)*256);
  int i;
  for (i=0;i<256;i++)
    img[i]=i;
  printf("finish\n");
  char fname[100];
  sprintf(fname,"screentest.fits");
  saveimagef(fname,img,16,16);
  free(img);
  */
  // test mkcmb
  /*  char fname[100];
  float *cl_l,*cl_p,*cl_py2,yp1,ypn;
  int cllength;
  long seed;
  struct simmap imap;
  initializemap(&imap);
  seed = -19;
  clinput(&cl_l,&cl_p,&cl_py2,&yp1,&ypn,&cllength);
  mkcmb(cllength,&seed,&imap,imap.img,cl_l,cl_p,cl_py2,yp1,ypn);
  sprintf(fname,"cmb.fits");
  savemap(&imap,fname);
  free(cl_l);
  free(cl_p);
  free(cl_py2);
  cleanmap(&imap);
  */
  float *img;
  long npix,i;
  char fname[200];
  npix=NDIM*NDIM;
  img = (float*)malloc(npix*sizeof(npix));
  for (i=0;i<npix;i++)
    img[i]=0;//1e-5;
/*  for (i=15;i<30;i++)
  {
    img[i+NDIM]=1e-4;
    img[i+2*NDIM]=1e-4;
    img[i+3*NDIM]=1e-4;
  }
*/  sprintf(fname,"sznull.fits");
  saveimagef(fname,img,NDIM,NDIM);
  return 0;
}
