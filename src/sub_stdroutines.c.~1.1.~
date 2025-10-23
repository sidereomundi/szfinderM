// $Id: sub_stdroutines.c,v 1.1 2012/06/04 14:03:20 jiayi Exp $
// changes: deleted unused variables

// write image as .imh or .fits 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "longnam.h"
#include "fitsio.h"

void saveimagef(char *name, float *data, int xdim, int ydim)
{
	int	i;
        fitsfile *fptr;
        int     status=0,ifchk_err();
        long    fpixel,nelements;
        int     bitpix=ULONG_IMG;
        long    numaxis=2,naxes[2];
        void    printerror(int);
        float   datamin=1.0e+30,datamax=-1.0e+30,bscale,bzero;

	// check to see if input image is fits image 
        if ((!strncmp("fits",&(name[strlen(name)-4]),4))
		 || (!strncmp("FITS",&(name[strlen(name)-4]),4))) {
	  // remove file if it exists 
          remove(name);
          // create the output file
          if (fits_create_file(&fptr,name,&status))
               printerror(status);
	  // set the axis length
          naxes[0]=xdim; naxes[1]=ydim; numaxis=2;
	  // set image to 32 bits
          bitpix=ULONG_IMG;
          if (fits_create_img(fptr,bitpix,numaxis,naxes,&status))
               printerror(status);
          fpixel=1;
          nelements=naxes[0]*naxes[1];
          // determine data range and bscale and bzero
          for (i=0;i<nelements;i++) {
            if (data[i]<datamin) datamin=data[i];
            if (data[i]>datamax) datamax=data[i];
          }
	  // this uses full dynamic range  (or am I wasting a factor of two?)
          bscale=(datamax-datamin)/2.0e+9;
	  // this sets minimum value to FITS min val
          bzero=datamin;
	  // save these new BZERO's and BSCALE's 
          if(fits_update_key(fptr,TFLOAT,"BSCALE",&bscale,
            "default scaling factor",&status))
            printerror(status);
          if(fits_update_key(fptr,TFLOAT,"BZERO",&bzero,
            "offset data range to that of unsigned short",&status))
            printerror(status);
          printf("  min %10.3e max %10.3e \n",datamin,datamax);
          // write the array of unsigned integers to the FITS file
          if (fits_write_img(fptr,TFLOAT,fpixel,nelements,data,&status) )
              printerror(status);
          // close the image
          if (fits_close_file(fptr,&status)) printerror(status);
	}
	else {
          printf("Only able to write FITS images, current name: %s\n",name);
          exit(1);
        }
}

void printerror(int status)
{
        fits_report_error(stderr,status);
        exit(1);
}

/* August 5th-- removed IRAF imh capability */
void readimage(char *name, float **dataptr,int *xdim, int *ydim)
{
        fitsfile *fptr;
        int     status=0, anynull,nfound;
        long    naxes[2],fpixel,npixels;
        float   *data,nullval;
        void    printerror(int);

	/* is image an IRAF imh file */
        if ((!strncmp("imh",&(name[strlen(name)-3]),3))
	    || (!strncmp("IMH",&(name[strlen(name)-3]),3))) {
	  printf("  IRAF images no longer supported\n");
	  exit(0);
	/*
          mode=IFRD_ONLY;
          imopen_(name,&mode,&im,&er,strlen(name));ifchk_err(er);
          imgsiz_(&im,axlen,&naxis,&dtype,&er);ifchk_err(er);
          *xdim=axlen[0];*ydim=axlen[1];
	  data=(float *)calloc((*xdim)*(*ydim),sizeof(float));
          for (i=1;i<=axlen[1];i++) {
            imgl2r_(&im,&(data[(i-1)*axlen[1]]),&i,&er);ifchk_err(er);
          }
          imclos_(&im,&er);ifchk_err(er);
	*/
	}
	/* is image a FITS file */
        else if ((!strncmp("fits",&(name[strlen(name)-4]),4))
		 || (!strncmp("fits",&(name[strlen(name)-4]),4))) {
          /* open the file */
          if (fits_open_file(&fptr,name,READONLY,&status)) printerror(status);
          /* read the NAXIS1 and NAXIS2 keyword to get image size */
          if (fits_read_keys_lng(fptr,"NAXIS",1,2,naxes,&nfound,&status))
            printerror(status);
	  *xdim=naxes[0];*ydim=naxes[1];
	  data=(float *)calloc((*xdim)*(*ydim),sizeof(float));
          /* prepare to read the image */
          npixels  = naxes[0]*naxes[1];
          fpixel   = 1; nullval  = 0;
          /* read the image */
          if (fits_read_img(fptr,TFLOAT,fpixel,npixels,&nullval,
            data,&anynull,&status)) printerror(status);
          /* close the image */
          if (fits_close_file(fptr,&status)) printerror(status);
	}
	else {
          printf("Only able to read FITS images\n");
          exit(0);
        }
	*dataptr=data;
}

// generate the filter size [arcmin]
double thetasize(int i)
{
	return (0.25+(float)i*0.25);
}
