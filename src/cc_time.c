/*
 *  Do correlation in time domain
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
#include <math.h>
#include "sacio.h"
#include "const.h"


void usage(void);

void usage() {
    fprintf(stderr, "Do correlation in time domain.                     \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Usage:                                             \n");
    fprintf(stderr, "  cc_time -Ttmark/ts/te [-h] sacfile1 sacfile2   \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Options:                                           \n");
    fprintf(stderr, "  -T: tmark/begin time (sec)/end timew (sec)       \n");
    fprintf(stderr, "  -h: show usage                                   \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Examples:                                          \n");
    fprintf(stderr, "  cc_time -T1/-3/4 seis1.sac seis2.sac             \n");
}

int main(int argc, char *argv[])
{
    int c;
    int error;
    int tmark;
    float t0, t1;

    float *data1, *data2;
    SACHEAD hd1, hd2;

    int i;
    double sum, sum1, sum2;

    error = 0;
    while ((c = getopt(argc, argv, "T:h")) != -1) {
        switch(c) {
            case 'T':
                if (sscanf(optarg, "%d/%f/%f", &tmark, &t0, &t1) != 3) error++;
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

    /* do correlation in time domain */
    sum  = 0.0;
    sum1 = 0.0;
    sum2 = 0.0;
    for (i=0; i<hd1.npts && i<hd2.npts; i++) {
        sum  += data1[i] * data2[i];
        sum1 += data1[i] * data1[i];
        sum2 += data2[i] * data2[i];
    }
    sum = sum / sqrt(sum1 * sum2);
    fprintf(stderr, "%s %s %g\n", argv[optind], argv[optind+1], sum);

    return 0;
}

