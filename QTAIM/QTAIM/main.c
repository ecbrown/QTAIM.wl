#include "mathlink.h"
#include <stdio.h>
#include <stdlib.h>

extern void read_wfn_init_(
int *filenamelength, 
char* filename, 
int nmonpnn[]
);

extern void read_wfn_(
int *filenamelength,
char* filename,
int *nmo,
int *np,
int *nn,
double xn[],
double zn[],
int centre[],
int type[],
double a[],
double o[],
double oe[],
double c[],
double ev[]
);

void readwfn(
int filenamelength, 
char filename[]
) {

  int nmonpnn[3];
  int nmo;
  int np;
  int nn;

  int len = filenamelength;

  read_wfn_init_(&len, filename, nmonpnn);

  nmo=nmonpnn[0];
  np=nmonpnn[1];
  nn=nmonpnn[2];

  //int filename_length;
  //char filename;
  //nmo;
  //np;
  //nn;
  double *xn=malloc(nn*3*sizeof(double));
  double *zn=malloc(nn*sizeof(double));
  int *centre=malloc(np*sizeof(int));
  int *type=malloc(np*sizeof(int));
  double *a=malloc(np*sizeof(double));
  double *o=malloc(nmo*sizeof(double));
  double *oe=malloc(nmo*sizeof(double));
  double *c=malloc(nmo*np*sizeof(double));
  double *ev=malloc(2*sizeof(double));

  read_wfn_(&len,filename,&nmo,&np,&nn,xn,zn,centre,type,a,o,oe,c,ev);


  MLPutFunction(stdlink, "List", 12);
   MLPutInteger32(stdlink, nmo);
   MLPutInteger32(stdlink, np);
   MLPutInteger32(stdlink, nn);
   MLPutReal64List(stdlink, xn, nn*3);
   MLPutReal64List(stdlink, zn, nn);
   MLPutInteger32List(stdlink, centre, np);
   MLPutInteger32List(stdlink, type, np);
   MLPutReal64List(stdlink, a, np);
   MLPutReal64List(stdlink, o, nmo);
   MLPutReal64List(stdlink, oe, nmo);   
   MLPutReal64List(stdlink, c, nmo*np);
   MLPutReal64List(stdlink, ev, 2);

  free(xn);
  free(zn);
  free(centre);
  free(type);
  free(a);
  free(o);
  free(oe);
  free(c);
  free(ev);

}

int main(int argc, char *argv[]) {
    return MLMain(argc, argv);
}

