/*********************************************************************
*	sac_wfcc.c:
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
*		05/07/03	Intepolate the output time shift to be below the samping point
*				by zpeng
*		08/18/20	change old sacio.c to new sacio.c by Jiayuan Yao
*
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
  float     il,ic,ih,maxl,maxc,maxh;

  double A[3][3], B[3];

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
        case 'D':
          sscanf(&argv[i][2],"%f/%f/%f",&tBefore,&tAfter,&tshift);
          break;
        case 'W':
          overWrite = 1;
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
       error = 1;
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

    /*update with more accurate determination of the time difference */
    maxc = crl[shift];

/*      printf("%d %e %e %d\n",shift,maxc,norm,max_shift); */

    if (norm != (mm/2-1)) {
      maxh = crl[shift+1];
    }
    else {
      maxh = crl[0];
    }
    if (norm !=0) {
      maxl = crl[shift-1];
    }
    else {
      maxl = crl[mm/2-1];
    }
    free(crl);

    il = -1.0; ic = 0.0; ih = 1.0;

    /* add on zpeng */
/*      printf("%e %e %e\n",maxl,maxc,maxh); */

    A[0][0]=1.0; A[0][1]=il; A[0][2]=il*il; B[0]=maxl;
    A[1][0]=1.0; A[1][1]=ic; A[1][2]=ic*ic; B[1]=maxc;
    A[2][0]=1.0; A[2][1]=ih; A[2][2]=ih*ih; B[2]=maxh;

    gauss(A,B,3,3,1.0e-6,&i,TRUE);
    ic = -B[1]/(2.0*B[2]);
    max = B[0] + B[1]*ic + B[2]*ic*ic;

    if( ic<(-1.0) || (ic>1.0) ) {
      fprintf(stderr,"warning: fractional correction ic = %f more than 1 sample\n",ic);
      exit(-1);
    }


/*      printf ("%f %e %d\n",ic,max,shift); */
    ic += shift;

    norm = max/modMaster/modOther;
    ic -= max_shift/2;
    arr -= ic*dt;

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

gauss(a,vec,n,nstore,test,ierror,itriag)
double *a, vec[], test;
int n, nstore, *ierror, itriag;
{

/* subroutine gauss, by william menke */
/* july 1978 (modified feb 1983, nov 85) */

/* a subroutine to solve a system of n linear equations in n unknowns*/
/* where n doesn't exceed MAX_TABLE_COLS */
/* gaussian reduction with partial pivoting is used */
/* 	a		(sent, destroyed)	n by n matrix		*/
/*	vec		(sent, overwritten)	n vector, replaced w/ solution*/
/* 	nstore		(sent)			dimension of a	*/
/*	test		(sent)			div by zero check number*/
/*	ierror		(returned)		zero on no error*/
/*	itriag		(sent)			matrix triangularized only*/
/*						 on TRUE useful when solving*/
/*						 multiple systems with same a */	static int isub[100], l1;
  int line[100], iet, ieb, i, j, k, l, j2;
  double big, testa, b, sum;


  iet=0;	/* initial error flags, one for triagularization*/
  ieb=0;	/* one for backsolving */

/* triangularize the matrix a*/
/* replacing the zero elements of the triangularized matrix */
/* with the coefficients needed to transform the vector vec */

  if (itriag) {	/* triangularize matrix */

     for( j=0; j<n; j++ ) {      /*line is an array of flags*/
      line[j]=0;
      /* elements of a are not moved during pivoting*/
      /* line=0 flags unused lines */
      }    /*end for j*/

    for( j=0; j<n-1; j++ ) {
      /*  triangularize matrix by partial pivoting */
                       big = 0.0; /* find biggest element in j-th column*/
          /* of unused portion of matrix*/
           for( l1=0; l1<n; l1++ ) {
                               if( line[l1]==0 ) {
                                       testa=(double) fabs(
            (double) (*(a+l1*nstore+j)) );
                                       if (testa>big) {
            i=l1;
            big=testa;
            } /*end if*/
          } /*end if*/
        } /*end for l1*/
                       if( big<=test) {   /* test for div by 0 */
                               iet=1;
                               } /*end if*/

                       line[i]=1;  /* selected unused line becomes used line */
                       isub[j]=i;  /* isub points to j-th row of tri. matrix */

                       sum=1.0/(*(a+i*nstore+j));
        /*reduce matrix towards triangle */
           for( k=0; k<n; k++ ) {
        if( line[k]==0 ) {
          b=(*(a+k*nstore+j))*sum;
                 for( l=j+1; l<n; l++ ) {
                                               *(a+k*nstore+l)=
              (*(a+k*nstore+l))
              -b*(*(a+i*nstore+l));
                 } /*end for l*/
                                       *(a+k*nstore+j)=b;
          } /*end if*/
        } /*end for k*/
      } /*end for j*/

               for( j=0; j<n; j++ ) {
      /*find last unused row and set its pointer*/
      /*  this row contians the apex of the triangle*/
      if( line[j]==0) {
        l1=j;   /*apex of triangle*/
        isub[n-1]=j;
        break;
        } /*end if*/
      } /*end for j*/

    } /*end if itriag true*/

  /*start backsolving*/

  for( i=0; i<n; i++ ) {	/* invert pointers. line(i) now gives*/
        /* row no in triang matrix of i-th row*/
        /* of actual matrix */
    line[isub[i]] = i;
    } /*end for i*/

  for( j=0; j<n-1; j++) { /*transform the vector to match triang. matrix*/               b=vec[isub[j]];
               for( k=0; k<n; k++ ) {
                      if (line[k]>j) {	/* skip elements outside of triangle*/
                                vec[k]=vec[k]-(*(a+k*nstore+j))*b;
        } /*end if*/
      } /*end for k*/
    } /*end for j*/

      b = *(a+l1*nstore+(n-1));   /*apex of triangle*/
      if( ((double)fabs( (double) b))<=test) {
    /*check for div by zero in backsolving*/
    ieb=2;
    } /*end if*/
      vec[isub[n-1]]=vec[isub[n-1]]/b;

      for( j=n-2; j>=0; j-- ) {	/* backsolve rest of triangle*/
    sum=vec[isub[j]];
    for( j2=j+1; j2<n; j2++ ) {
      sum = sum - vec[isub[j2]] * (*(a+isub[j]*nstore+j2));
      } /*end for j2*/
      b = *(a+isub[j]*nstore+j);
               if( ((double)fabs((double)b))<=test) {
      /* test for div by 0 in backsolving */
      ieb=2;
      } /*end if*/
    vec[isub[j]]=sum/b;   /*solution returned in vec*/
    } /*end for j*/

/*put the solution vector into the proper order*/

      for( i=0; i<n; i++ ) {    /* reorder solution */
    for( k=i; k<n; k++ ) {  /* search for i-th solution element */
      if( line[k]==i ) {
        j=k;
        break;
        } /*end if*/
      } /*end for k*/
               b = vec[j];       /* swap solution and pointer elements*/
               vec[j] = vec[i];
               vec[i] = b;
               line[j] = line[i];
    } /*end for i*/

    *ierror = iet + ieb;   /* set final error flag*/
}

