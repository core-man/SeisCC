/*
 *  Do correlation in frequency domain using SAC's libsac.a, i.e., crscor
 *
 *  Author: Jiayuan Yao @ NTU
 *
 *  Revisions:
 *    2020-08-18  Jiayuan Yao   revised from SAC's example: correlatec.c
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <malloc.h>
#include <math.h>
#include "sacio.h"
#include "sac.h"
#include "const.h"

void usage();

void usage()
{
    fprintf(stderr, "Do correlation in frequency domain calling crscor  \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Usage:                                             \n");
    fprintf(stderr, "  correlatec -Ttmark/ts/te [-Occf] [-Wtaper] [-Acczero]  \n");
    fprintf(stderr, "             [-h] sacfile1 sacfile2                \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Options:                                           \n");
    fprintf(stderr, "  -T: tmark/begin time (sec)/end time (sec)        \n");
    fprintf(stderr, "  -O: cross-correlation function file              \n");
    fprintf(stderr, "  -W: taper (0: NO (default); 1: hanning; 2: cos)            \n");
    fprintf(stderr, "  -A: only output cross-correlation at zero lag time (0: YES; 1: NO (default)) \n");
    fprintf(stderr, "  -h: show usage                                   \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Examples:                                          \n");
    fprintf(stderr, "  correlatec -T1/-3/4 seis1.sac seis2.sac            \n");
    fprintf(stderr, "  correlatec -T1/-3/4 -A0 seis1.sac seis2.sac        \n");
    fprintf(stderr, "  correlatec -T1/-3/4 -W1 seis1.sac seis2.sac        \n");
    fprintf(stderr, "  correlatec -T1/-3/4 -Occf12.sac seis1.sac seis2.sac\n");
}


int main(int argc, char *argv[])
{
    int c,err;
    int tmark;
    float t0, t1;
    int taper=0, cczero=1;
    char outfile[MAX_FNAME];

    err = 0;
    while ((c = getopt(argc, argv, "T:W:O:A:h")) != -1) {
        switch(c) {
            case 'T':
                if (sscanf(optarg, "%d/%f/%f", &tmark, &t0, &t1) != 3) err++;
                break;
            case 'W':
                if (sscanf(optarg, "%d", &taper) != 1) err++;
                break;
            case 'O':
                if (sscanf(optarg, "%s", outfile) != 1) err++;
                break;
            case 'A':
                if (sscanf(optarg, "%d", &cczero) != 1) err++;
                break;
            case 'h':
                usage();
                return -1;
            default:
                usage();
                return -1;
        }
    }

    if (argc-optind != 2 || err) {
        usage();
        return -1;
    }


    /* Local variables */
    int i;
    int nlen, nlen1, nlen2, max;
    /*int nerr;*/

    float beg1, beg2, delta;
    /*float end;*/
    /*char *kname;*/

    float *yarray1, *yarray2, *ytmp;
    double *dtap;
    SACHEAD hd1, hd2, hd;
    float *out;
    float max_value, max_time;

    int nwin, wlen, nfft;
    /*int leven, when; */

    double pw1, pw2, scale;

    char error[ERROR_MAX];

    /*max = MAX;*/


    /* read data and check delta */
    yarray1 = read_sac_pdw(argv[optind],   &hd1, tmark, t0, t1);
    yarray2 = read_sac_pdw(argv[optind+1], &hd2, tmark, t0, t1);
    if (fabs(hd1.delta - hd2.delta) > EPS) {
        fprintf(stderr, "delta is not equal: %f %f.\n", hd1.delta, hd2.delta);
        return -1;
    }
    nlen1 = hd1.npts;
    beg1  = hd1.b;
    delta = hd1.delta;
    nlen2 = hd2.npts;
    beg2  = hd2.b;

    /* Read in the first file  */
    /*rsac1(kname, yarray1, &nlen1, &beg1, &delta, &max, &nerr, SAC_STRING_LENGTH);
    rsac1(kname, yarray2, &nlen2, &beg2, &delta, &max, &nerr, SAC_STRING_LENGTH);*/

    nlen = nlen1;
    /**
     *  If the signals are not the same length, then find the longest
     *  signal, make both signals that length by filling the remainder
     *  with zeros (pad at the end) and then run them through crscor
     *  This should be fixed in upcoming releases and the introduction
     *  of a "correlate" function so you do not need to handle the
     *  signal length, padding, and window lengths, ...
     *
     */


    /* set memory for taper */
    if ((dtap = malloc(nlen * sizeof(*dtap))) == NULL) {
        fprintf(stderr, "Allocation failed for taper.\n");
        exit(-1);
    }

    /* taper data */
    if (taper == 1) {
        taper_hanning(dtap, nlen, 0.05);  /* happing */
    } else if (taper == 2) {
        taper_cos(dtap, nlen, 0.05);      /* cos     */
    }
    if (taper == 1 || taper ==2) {
        for (i = 0; i < nlen; i++) {
            yarray1[i] = yarray1[i] * dtap[i];
            yarray2[i] = yarray2[i] * dtap[i];
        }
    }


    /* Allocate space for the correlation of yarray1 and yarray2 */
    max = next2((2 * nlen) - 1) * 2;
    if((out = (float *) malloc(sizeof(float) * max)) == NULL) {
      fprintf(stderr, "Error allocating memory for correlation\n");
      exit(-1);
    }
    /* Allocate space for ytmp */
    if ((ytmp = (float *) malloc(sizeof(*ytmp) * 4 * nlen)) == NULL) {
      fprintf(stderr, "Allocation failed for ytmp\n");
      exit(-1);
    }

    /* normalization factor */
    for (pw1 = 0.0, i = 0; i < nlen1; i++)
        pw1 += yarray1[i] * yarray1[i];
    for (pw2 = 0.0, i = 0; i < nlen2; i++)
        pw2 += yarray2[i] * yarray2[i];
    scale = 1.0 / sqrt(pw1 * pw2);

    /* Set up values for the cross correlation */
    nwin = 1;
    wlen = nlen;
    nfft = 0;

    /*     Call crscor ( Cross Correlation )
     *        - yarray1 - First  Input array to correlate
     *        - yarray2 - Second Input array to correlate
     *        - nlen    - Number of points in yarray and yarray2
     *        - nwin    - Windows to use in the correlation
     *        - wlen    - Length of the windows
     *        - type    - Type of Window (SAC_RECTANGLE)
     *        - out     - output sequence (2*wlen-1)
     *        - nfft    - Length of the output sequence
     *        - error   - Error Message
     *        - err_len - Length of Error Message (on input)
     */
    crscor(yarray1, yarray2, nlen,
           nwin, wlen, SAC_RECTANGLE,
           ytmp, &nfft, error, ERROR_MAX);

    for(i = 0; i < max; i++) {
      out[i] = 0;
    }

    /*
     *     out[0 : nlen1 - 2 ] <-- ytmp[ nfft - nlen1 + 1 : nfft -1 ]
     *     out[nlen1 - 1 : nlen1 + nlen2 - 2 ] <-- ytmp[ 0 : nlen2-1 ]
     */
    for(i = 0; i <= nlen1 - 2; i++) {
      out[i] = ytmp[nfft - nlen1 + i + 1] * scale;
    }
    for(i = 0; i <= nlen2 - 1; i++) {
      out[nlen1 + i - 1] = ytmp[i] * scale;
    }


    nfft      = nlen1 + nlen2 - 1;
    /*xarray[0] = 0;
    leven     = TRUE;*/
    beg1      = -delta * (nlen1 - 1) + (beg2 - beg1);
    /*end       = beg1 + delta * (nfft - 1);*/

    /*
    setnhv ( "npts",   &nfft,   &nerr, SAC_STRING_LENGTH);
    setfhv ( "delta",  &delta,  &nerr, SAC_STRING_LENGTH);
    setlhv ( "leven",  &leven,  &nerr, SAC_STRING_LENGTH);
    setfhv ( "b",      &beg1,   &nerr, SAC_STRING_LENGTH);
    setfhv ( "e",      &end,    &nerr, SAC_STRING_LENGTH);
    setihv ( "iftype", "itime", &nerr, SAC_STRING_LENGTH, SAC_STRING_LENGTH);
    when = SAC_NUMBER_UNDEFINED;
    setnhv ( "nzyear", &when,   &nerr, SAC_STRING_LENGTH);
    setnhv ( "nzhour", &when,   &nerr, SAC_STRING_LENGTH);
    *
    */

    hd = new_sac_head(hd1.delta, nfft, beg1);

    /* Find the maximum value and time of the correlation function */
    //max_time  = 0;
    max_time  = (1 - nlen) * delta;
    max_value = out[0];
    for(i = 1; i < nfft; i++) {
      if(fabs(out[i]) > max_value) {
        max_value = fabs(out[i]);
        /* Negative shifts are at the end of the correlation sequence */
        /*
        if(i > nfft/2) {
          max_time  = (i - nfft) * delta;
        } else {
          max_time  = i * delta;
        }*/
        max_time = (i + 1 - nlen) * delta;
      }
    }

    /*
      setfhv( "user0", &max_time,  &nerr, SAC_STRING_LENGTH);
      setfhv( "user1", &max_value, &nerr, SAC_STRING_LENGTH);
    */

    /* output cross-correlation or auto-correlation coefficient and time shift */
    if (cczero != 0) {
        fprintf(stdout, "%s %s %f %f\n", argv[optind], argv[optind+1], max_value, (hd2.b-hd1.b)+max_time);
    } else {
        fprintf(stdout, "%s %s %f\n", argv[optind], argv[optind+1], out[nlen-1]);
    }

    /*   Write out the correlation function   */
    if (strlen(outfile) != 0)
        write_sac(outfile, hd, out);


    /* free memory */
    free(yarray1); free(yarray2); free(ytmp);
    free(dtap);
    free(out);

    return 0;
}

