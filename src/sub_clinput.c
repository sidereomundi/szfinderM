//$Id: sub_clinput.c,v 1.4 2012/07/10 13:08:51 jiayi Exp $
// read in the cl.dat ( From Amy Lien )

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parameters.h"
#include "nr.h"

#define clswitch 1  // switch for old/new cl file (0/1)

void clinput(float **cl_l,float **cl_p,float **cl_py2,float *yp1,float *ypn,int *cllength)
{
  FILE	*pip,*inp;
  char	trash[200],command[100];
  int	i;

//  FILE *ftest;
//  ftest=fopen("cltest.dat","w");
#if clswitch == 0
  printf(" read in old version of cl.dat \n");
#else
  printf(" read in new/CAMB version of cl.dat \n");
#endif

  /* read in CMB power spectrum */
  /* first calculate file size*/
  sprintf(command,"wc %s",clfile);
  pip=popen(command,"r");
  fgets(trash,200,pip);
  pclose(pip);
  sscanf(trash,"%d",cllength);
  inp=fopen(clfile,"r");
  if(inp==NULL)
    {
      printf("# fail to open %s\n",clfile);
      exit(1);
    }
  /* then allocate array space */
  *cl_l=(float *)calloc(*cllength,sizeof(float));
  *cl_p=(float *)calloc(*cllength,sizeof(float));
  *cl_py2=(float *)calloc(*cllength,sizeof(float));
  /* then actually read in file */
  i=1;
  while (fgets(trash,200,inp)!=NULL) {
    // read in the CAMB output powerspectrum
    // the difference is demonstrated in checkcl.pro
    sscanf(trash,"%f %f",(*cl_l)+i,(*cl_p)+i);
    /* store the rms in dT in cl_p vector */
    (*cl_p)[i] *= (2.0*M_PI/((*cl_l)[i]*((*cl_l)[i]+1.0)));
#if clswitch == 0
    (*cl_p)[i]=T_CMB*sqrt((*cl_p)[i]); // (old version)
#else
    (*cl_p)[i] *= 1.e-12;  // in uK (CAMB result)
    (*cl_p)[i]=sqrt((*cl_p)[i]);
#endif
//    fprintf(ftest,"%g %g\n",(*cl_l)[i],(*cl_p)[i]);
    i++;
  }
  fclose(inp);
//  fclose(ftest);
  
  printf("## %d from %d lines from cl.dat\n",i-1,*cllength); fflush(stdout);
  /* create spline of c_l's for future interpolation */
  /* calculate derivatives at first and last elements */
  *yp1=((*cl_p)[2]-(*cl_p)[1])/((*cl_l)[2]-(*cl_l)[1]);
  *ypn=((*cl_p)[*cllength]-(*cl_p)[*cllength-1])/((*cl_l)[*cllength]-(*cl_l)[*cllength-1]);
  spline(*cl_l,*cl_p,*cllength,*yp1,*ypn,*cl_py2);
}
