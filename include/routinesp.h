//$Id: routinesp.h,v 1.2 2013/02/08 10:34:56 jiayi Exp $
// routines for for c++
#ifndef _ROUTINESP_H_
#define _ROUTINESP_H_
#include <vector>
using namespace std;
struct geometry;
struct peak;
struct simcluster;

// sub_cluster.cpp
void loadgeometry(char *, struct geometry *);  // load wmap.geometry.dat to find the redshift for clusters
float findredshift(int , struct geometry *);   // find the cluster redshift in geometry file
void inputcluster(char *, vector<struct simcluster> *);

// sub_peak.cpp
void peaklist(int,vector<struct peak> *,float *);  // find peak list
void peakcombine(vector<struct peak> *);  // remove duplicate peaks
void peakoutput(vector<struct peak> *,char*);
void peakinput(vector<struct peak> *,char*);
#endif
