/*
 *  Do correlation in frequency domain
 *
 *  Author: Jiayuan Yao @ NTU
 *
 *  Revisions:
 *    2020-08-17  Jiayuan Yao  Initial Coding
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "complex.h"
#include "fftw3.h"
#include "sacio.h"
#include "const.h"


typedef struct _ccvalue{
    int     imax;
    double ccmax;
    double scale;
} ccvalue;

void usage(void);
ccvalue cc_freq(int nd, int nd1, float *data, float *data1, double *cc, int taper);


int main(int argc, char *argv[])
{
    int c,error;
    int tmark;
    float t0, t1;
    int taper=0, autoc=1;
    char outfile[MAX_FNAME];

    float *data1, *data2;
    SACHEAD hd1, hd2, hd;

    int i, cc_npts;
    float b, *cc_norm;
    double *cc;
    ccvalue ccv;

    error = 0;
    while ((c = getopt(argc, argv, "T:W:O:A:h")) != -1) {
        switch(c) {
            case 'T':
                if (sscanf(optarg, "%d/%f/%f", &tmark, &t0, &t1) != 3) error++;
                break;
            case 'W':
                if (sscanf(optarg, "%d", &taper) != 1) error++;
                break;
            case 'O':
                if (sscanf(optarg, "%s", outfile) != 1) error++;
                break;
            case 'A':
                if (sscanf(optarg, "%d", &autoc) != 1) error++;
                break;
            case 'h':
                usage();
                return -1;
            default:
                usage();
                return -1;
        }
    }

    if (argc-optind != 2 || error) {
        usage();
        return -1;
    }


    /* read data and check delta */
    data1 = read_sac_pdw(argv[optind], &hd1, tmark, t0, t1);
    data2 = read_sac_pdw(argv[optind+1], &hd2, tmark, t0, t1);
    if (fabs(hd1.delta - hd2.delta) > EPS) {
        fprintf(stderr, "delta is not equal: %f %f.\n", hd1.delta, hd2.delta);
        return -1;
    }

    cc_npts = hd1.npts + hd2.npts - 1;

    /* set memory for cross-correlation function */
    if ((cc = malloc(cc_npts * sizeof(*cc))) == NULL) {
        fprintf(stderr, "Allocation failed for cross-correlation function\n");
        return -1;
    }

    /* set memory for normalized cross-correlation function */
    if ((cc_norm = malloc(cc_npts * sizeof(*cc_norm))) == NULL) {
        fprintf(stderr, "Allocation failed for normalized cross-correlation function\n");
        return -1;
    }

    /* do cross-correlate in frequency domain  */
    ccv = cc_freq(hd1.npts, hd2.npts, data1, data2, cc, taper);

    /* output cross-correlation or auto-correlation coefficient and time shift */
    if (autoc != 0) {
        fprintf(stdout, "%s %s %f %f\n", argv[optind], argv[optind+1],
                ccv.ccmax*ccv.scale, ((hd2.b - hd1.b) - ccv.imax*hd1.delta));
    } else {
        fprintf(stdout, "%s %s %f\n", argv[optind], argv[optind+1],
                cc[hd1.npts-1]*ccv.scale);
    }

    /* get cross-correlation function */
    if (strlen(outfile) != 0) {
        for(i = 0; i < cc_npts; i++)
             cc_norm[i] = cc[i] * ccv.scale;

        //hd.b = (-hd1.npts + 1) * hd1.delta + (hd2.ts - hd1.ts);
        b = (-hd1.npts + 1)*hd1.delta + (hd2.b - hd1.b);
        hd = new_sac_head(hd1.delta, cc_npts, b);
        //strcpy(hd.kstnm,  hd1.kstnm);
        //strcpy(hd.knetwk, hd1.knetwk);
        //strncpy(hd.kstnm,  hd1.kstnm,  8);
        //strncpy(hd.knetwk, hd1.knetwk, 8);

        write_sac(outfile, hd, cc_norm);
    }

    free(data1);
    free(data2);
    free(cc);
    free(cc_norm);


    return 0;
}



/*
 * cross-correlation in frequency domain
 */
ccvalue cc_freq(int nd1, int nd2, float *data, float *data1, double *cc, int taper)
{
    int k, num;
    int nw;
    double *d1, *d2, *dtap;
    double pw1, pw2;

    int imax;
    double ccmax, scale, cc0;
    ccvalue ccv;

    int nfft;
    double inv_nfft;
    double complex *fft_in, *fft_out, *fft_out1;
    fftw_plan pf, pb;


    /* set time length and nfft */
    if (nd1 >= nd2) {
        for (k=2*nd1,  nfft=1; nfft < k; nfft += nfft);
        nw = nd1;
    }
    else {
        for (k=2*nd2, nfft=1; nfft < k; nfft += nfft);
        nw = nd2;
    }
    inv_nfft = 1.0 / (double) nfft;

    /* set memory for data (first SAC, second SAC, taper */
    if ((d1 = malloc(3 * nw * sizeof(*d1))) == NULL) {
        fprintf(stderr, "Allocation failed for data.\n");
        exit(-1);
    }
    d2   = d1 + nw;   /* second SAC data */
    dtap = d2 + nw;   /* taper           */

    /* taper data */
    if (taper == 1) {
        taper_hanning(dtap, nw, 0.05);  /* happing */
    } else if (taper == 2) {
        taper_cos(dtap, nw, 0.05);      /* cos     */
    }

    for (k = 0; k < nd1; k++) {
        if (taper == 1 || taper ==2) {
            d1[k] = data[k] * dtap[k];
        } else {
            d1[k] = data[k];
        }
    }
    for (k = 0; k < nd2; k++) {
        if (taper == 1 || taper ==2) {
            d2[k] = data1[k] * dtap[k];
        } else {
            d2[k] = data1[k];
        }
    }

    /* normalization factor */
    for (pw1 = 0.0, k = 0; k < nd1; k++)
        pw1 += d1[k] * d1[k];
    for (pw2 = 0.0, k = 0; k < nd2; k++)
        pw2 += d2[k] * d2[k];
    scale = 1.0 / sqrt(pw1 * pw2);


    /* do FFT */
    if ((fft_in = malloc(3 * nfft * sizeof(*fft_in))) == NULL) {
        fprintf(stderr, "Allocation failed for fft_in\n");
        exit(-1);
    }

    fft_out  = fft_in + nfft;
    fft_out1 = fft_out + nfft;
    pf = fftw_plan_dft_1d(nfft, fft_in, fft_out,  FFTW_FORWARD,  FFTW_ESTIMATE);
    pb = fftw_plan_dft_1d(nfft, fft_in, fft_out1, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (k = 0; k < nfft; k++) {
        if(k < nd1)
            fft_in[k] = d1[k];
        else
            fft_in[k] = 0.0;
    }
    fftw_execute_dft(pf, fft_in, fft_out);

    for (k = 0; k < nfft; k++) {
        if (k < nd2)
            fft_in[k] = d2[k];
        else
            fft_in[k] = 0.0;
    }
    fftw_execute_dft(pf, fft_in, fft_out1);

    for (k = 0; k < nfft; k++)
        fft_out[k] *= conj(fft_out1[k]) * inv_nfft;
    fftw_execute_dft(pb, fft_out, fft_in);
    fftw_destroy_plan(pf);
    fftw_destroy_plan(pb);

    /* get cross-corrleation function */
    imax  = 0;
    ccmax = creal(fft_in[0]);
    cc0   = fabs(ccmax);
    /* cc[nd1-1 : 0]  <---  fft_in[0 : nd1-1]
     * nd1-1 is zero lag time */
    for (num = 0, k = nd1-1; num < nd1; num++) {
        cc[k] = creal(fft_in[num]);
        if (fabs(cc[k]) > cc0) {
            imax  = num;
            ccmax = cc[k];
            cc0   = fabs(cc[k]);
        }
        k--;
    }
    /* cc[nd1 : nd1+nd2-2]  <---  fft_in[nd1+nd2-2 : nd1] */
    for (num = 1, k = nd1; num < nd2 ; num++) {
        cc[k] = creal(fft_in[nfft - num]);
        if (fabs(cc[k]) > cc0) {
            imax  = -num;
            ccmax = cc[k];
            cc0   = fabs(cc[k]);
        }
        k++;
    }
    free(d1); free(fft_in);

    /* get maximum cross-correlation coefficient */
    ccv.imax  = imax;
    ccv.ccmax = ccmax;
    ccv.scale = scale;

    return ccv;
}


void usage()
{
    fprintf(stderr, "Do correlation in frequency domain.                \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Usage:                                             \n");
    fprintf(stderr, "  cc_freq  -Ttmark/ts/te -Occf -Wtaper -Aautocorr  \n");
    fprintf(stderr, "           [-h] sacfile1 sacfile2                  \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Options:                                           \n");
    fprintf(stderr, "  -T: tmark/begin time (sec)/time window (sec)     \n");
    fprintf(stderr, "  -O: cross-correlation function file              \n");
    fprintf(stderr, "  -W: taper (0: NO; 1: hanning; 2: cos)            \n");
    fprintf(stderr, "  -A: output auto-correlation                      \n");
    fprintf(stderr, "  -h: show usage                                   \n");
}

