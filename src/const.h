/*******************************************************************************
    Name:     para.h

    Purpose:    some common variables

    Revisions:
        08/17/2020  Jiayuan Yao  Initial coding
*******************************************************************************/

#ifndef PARA_H
#define PARA_H

#define ERROR_MAX  256      /* maxmum error length        */
#define MAX_FNAME  256      /* maxmum file name length    */
#define EPS        10e-5    /* float number equality flag */

/* function prototype of taper functions */
int taper_cos(double *tpr, int n, double per);
int taper_hanning(double *tpr, int n, double per);

#endif /* para.h */
