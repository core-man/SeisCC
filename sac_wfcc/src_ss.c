/*********************************************************************
*	src_ss.c:
*		estimate teleseismic source function by stacking
*		all vertical components of an array. The code uses
*		cross-correlation between each trace and a master
*		trace to align all traces at their first arrivals.
*
*	Usage:
*		src_ss -Dlength_of_correlation_window/max_time_shift [-W] [-N]
*		then input name of file from stdin, in following format:
*			name arr
*		where arr is the approx. arrival time for the trace.
*		The first line is for the master trace.
*
*		The outputs are in the same format so it can be used
*		as input for iteration
*
*	Author:  Lupei Zhu
*
*	Revision History
*		June 1997	Initial coding
*		06/23/97	input component names from stdio with
*				a shift0.
*		09/04/97	change shift0 to arrival time
*		02/11/02	output cross-correlation value
*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sacio.h"
#include "Complex.h"

#define TAPER	0.2		/* portion of window being tapered */

int main(int argc, char **argv) {
  int 		i, nn, mm, mt8, t8, max_shift, error;
  int		shift, start, end, ntrace, overWrite, normalize;
  char		line[128],inf[64],outf[64];
  size_t	wndw_size;
  float		tshift;	/*max. time shift in sec.*/
  float		tBefore, tAfter, arr, max, min, modMaster, modOther;
  float 	norm, dt, *src, *master, *other, *trace, *crl;
  SACHEAD	hd_m, hd;
  void		taper(float *, int);

  error = 0;
  tBefore = -5;
  tAfter = 10;
  tshift = 1;
  overWrite = 0;
  normalize = 0;
  /* input parameters */
  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
       switch(argv[i][1]) {

       case 'W':
	 overWrite = 1;
         break;

       case 'D':
         sscanf(&argv[i][2],"%f/%f/%f",&tBefore,&tAfter,&tshift);
         break;

       case 'N':
	 overWrite = 1;
	 normalize = 1;
	 break;

       default:
         error = 1;
         break;

       }
    } else {
       error  = 1;
    }
  }

  if (argc == 1 || error) {
     fprintf(stderr,"usage: %s -Dt1/t2/max_shift [-W] [-N]\n",argv[0]);
     return -1;
  }

  /* input master trace */
  gets(line);
  sscanf(line,"%s %f",outf,&arr);
  if ( (src=read_sac(outf,&hd_m)) == NULL ) return -1;
  printf("%-6s %8.4f %7.4f\n",outf,arr,1.0);
  fflush(stdout);
  nn = hd_m.npts;
  dt = hd_m.delta;
  max_shift = 2*rint(tshift/dt);
  mm = rint((tAfter-tBefore)/dt);
  wndw_size = mm*sizeof(float);
  if ( (master=(float *)malloc(wndw_size)) == NULL ||
       (other=(float *)malloc(wndw_size)) == NULL ) {
    fprintf(stderr,"fail to allocation memory for src\n");
    return -1;
  }
  mt8 = rint((arr+tBefore-hd_m.b)/dt);
  if (mt8 < 0) {
    fprintf(stderr,"%s time before arr. is not long enough\n",outf);
    return -1;
  }
  memcpy(master, src+mt8, wndw_size);
  taper(master, mm);
  for(modMaster=0.,i=0;i<mm;i++) modMaster += master[i]*master[i];
  modMaster = sqrt(modMaster);

  ntrace = 0;
  for(i=0;i<nn;i++) src[i] = 0.;
  while (gets(line)) {

    sscanf(line,"%s %f",inf,&arr);

    if ( (trace=read_sac(inf,&hd)) == NULL ) continue;
    t8 = rint((arr+tBefore-hd.b)/dt);
    if (t8 < 0) {
      fprintf(stderr,"%s time before arr. is not long enough\n",inf);
      continue;
    }
    memcpy(other, trace+t8,wndw_size);
    taper(other, mm);
    for(modOther=0.,i=0;i<mm;i++) modOther += other[i]*other[i];
    modOther = sqrt(modOther);

    crl = crscrl(mm, master, other, max_shift);
    shift = 0;
    norm = 0;
    for(i=0;i<=max_shift;i++) {
      if (crl[i]>norm) {
	shift = i;
	norm = crl[i];
      }
    }
    free(crl);

    norm = norm/modMaster/modOther;
    shift -= max_shift/2;
    arr -= shift*dt;
    printf("%-6s %8.4f %7.4f\n",inf,arr,norm);
    fflush(stdout);
    shift -= (t8 - mt8);

    /* stacking */
    if (! overWrite ) continue;
    ntrace++;
    start = shift;		if (start<0) start = 0;
    end   = hd.npts+shift;	if (end>nn) end = nn;
    norm = 1.;
    if (normalize) {
      min = 1.e+32;
      max =-1.e+32;
      for (i=0; i<mm; i++) {
        if (other[i]>max) max=other[i];
        if (other[i]<min) min=other[i];
      }
      norm = 1./(max-min);
    }
    for (i=start; i<end; i++)
      src[i] += norm*trace[i-shift];
    free(trace);
  }

  if (ntrace<1 || ! overWrite) return 0;
  norm = 1./ntrace;
  for(i=0;i<nn;i++) src[i] *= norm;
  write_sac(outf, hd_m, src);
  free(src);
  return 0;

}

void	taper(float *aa, int n)
{
  int i, m;
  float	tt, pi1;
  m = TAPER*n;
  pi1 = 3.1415926/m;
  for (i=0; i<m; i++) {
    tt = 0.5*(1.-cos(i*pi1));
    aa[i] *= tt;
    aa[n-i-1] *= tt;
  }
}
