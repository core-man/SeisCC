
This code is downloaded from [Introduction to Seismic Analysis Code (SAC)](http://geophysics.eas.gatech.edu/people/zpeng/Teaching/SAC_Tutorial/) on Aug. 18 2020. Below is original `readme`.


This directory contains the following two c programs `src_ss` and `sac_wfcc`, which
can be used to compute the cross-correlation coefficients between pairs of
seismograms. The main difference between `src_ss` and `sac_wfcc` is that `src_ss`
output the time alignment and correlation coefficient at the time interval
of the sampling point, which `sac_wfcc` does a simple interpolation to obtain
subsampling accuracy of time alignment and correlation coefficient.

The `src_ss` code is originally written by Dr. Lupei Zhu at SLU. The subroutine
to get the interpolation gauss is written by William Menke at Lamont.

```bash
Usage:
$ saclst a f *.NC.CCO.EHZ.SAC |  sac_wfcc -D-0.2/0.8/0.5
# list the original P wave arrival, and perform waveform cross-correlations using 0.2 s before
# and 0.8 s after the P arrival. The output is the wf, the P arrival after time shift to best align
# with the first trace, and the correlation coefficient.
# the related example can be downloaded from: http://geophysics.eas.gatech.edu/people/zpeng/Teaching/SAC_Tutorial/
```

Last updated by zpeng@gatech.edu, Wed Feb 24 21:59:33 EST 2010

