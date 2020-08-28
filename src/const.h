/*******************************************************************************
    Name:     const.h

    Purpose:    some constant variables

    Revisions:
        08/17/2020  Jiayuan Yao  Initial coding
*******************************************************************************/

#ifndef CONST_H
#define CONST_H

#define ERROR_MAX  256      /* maxmum error length        */
#define MAX_FNAME  256      /* maxmum file name length    */
#define EPS        10e-5    /* float number equality flag */

/* function prototype of taper functions */
int taper_cos(double *tpr, int n, double per);
int taper_hanning(double *tpr, int n, double per);

#endif /* const.h */
