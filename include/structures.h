//$Id: structures.h,v 1.5 2012/06/11 13:57:55 jiayi Exp $
// Define the structures used for sz matched filter
#ifndef _STRUCTURES_H_
#define _STRUCTURES_H_

#define NZSLICE 33
// simulation map
struct simmap
{
  double pixsize;  // pixel size [arcmin]
  unsigned int nsidepix;
  long npix; // total # of pixel
  float *img;
};

// peak detected from mock sky
struct peak
{
  int xpix,ypix;
  float sn,coresize;
  // float redshift;
  // int matched;
};

// cluster from simulation catalog
struct simcluster
{
  float xpos,ypos;   // pixel position from the simulation
  int xpix,ypix;
  int ihal;
  double m;          // viral mass [unit??]
  double rad;        // Rrel from simulation catalog
  float redshift;
  int flag;          // the flag for cluster center is in the lightcone slice (=0), else =1.
};

struct geometry
{
  int isnap[NZSLICE];
  double zslice[NZSLICE];
};
#endif
