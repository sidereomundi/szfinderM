//$Id: sub_mkbolo.c,v 1.2 2012/06/13 09:25:05 jiayi Exp $
// Generate SZ+CMB Temperature based on CMB temperature

#include <math.h>
#include <stdio.h>
#include "nr.h"
#include "structures.h"
#include "parameters.h"
#include "routines.h"

// Convert temperature Map to intensity map
// input note:  imap  - compton y map
//              noise - empty float array
//              Tsig  - SPT whitenoise in [K/arcmin]
void mkbolo(imap, frequency, cmb, noise, Tsig, seed)
  struct simmap *imap;
  double frequency, Tsig;
  float *cmb, *noise;
  long *seed;
{
  int j,size;
  double x,expx,temperature;
  double hvk;//,khc,pbk;
  double ftemp;

  // initial setup
  size=imap->npix;
  frequency*=1.0e9;       // [GHz] to [Hz]
  Tsig*=1.0e-6/imap->pixsize;   // origin unit [K arcmin^-1]  to unit [K/pixel]
  
  hvk = H*frequency/KB;
//  khc = 2*Cube(KB)/Squ(H*C_CGS);
//  pbk = 2*H*Cube(frequency)/Squ(C_CGS);
  for (j=0;j<size;j++) {
    // calculate sz + cmb flux
    // Birkinshaw 99's paper (not quite right...)
    /*	    	temperature=cmb[j];
		x=hvk/temperature;
		expx=exp(x);
		flux=1.0+sz[j]*Quad(x)*expx/Squ(expx-1.0)*(x*(expx+1.0)/(expx-1.0)-4.0);
		ftemp=pbk/(expx-1.0);
		sz[j] = flux*ftemp;
    */	
    // Klaus's way of define y
    temperature=cmb[j];
    x=hvk/temperature;
    expx=exp(x);
    ftemp=x*(expx+1.0)/(expx-1.0)-4.0;
    ftemp=(imap->img[j]*ftemp+1)*temperature;
    imap->img[j]=ftemp;
    //	sz[j] = 0; // for sz-null image (pure_noise.fits)
    // calculate white noise temperature
    noise[j]=gasdev(seed)*Tsig;
  }
}
