/*
 *  Some taper functions
 *
 *  Author: Jiayuan Yao @ NTU
 *
 *  Revisions:
 *    2020-08-17  Jiayuan Yao  Initial Coding
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int taper_cos(double *tpr, int n, double per)
{
    int i, nptr;
    double omega, v;

    if (per > 0.5)
        return -1;

    for (i = 0; i < n; i++)
        tpr[i] = 1.0;

    nptr = (int) ((n + 1) * per);
    if (nptr < 2)
        nptr = 2;
    omega = 0.5 * M_PI / nptr;

    for (i = 0; i < nptr; i++) {
        v = sin(omega * i);
        tpr[i] *= v;
        tpr[n - 1 - i] *= v;
    }

    return 0;
}


int taper_hanning(double *tpr, int n, double per)
{
    int i, nptr;
    double omega, v, f0 = 0.5, f1 = 0.5;

    if (per > 0.5)
        return -1;

    for (i = 0; i < n; i++)
        tpr[i] = 1.0;

    nptr = (int) ((n + 1) * per);
    if (nptr < 2)
        nptr = 2;
    omega = M_PI / nptr;

    for (i = 0; i < nptr; i++) {
        v = f0 - f1 * cos(omega * i);
        tpr[i] *= v;
        tpr[n - 1 - i] *= v;
    }

    return 0;
}


