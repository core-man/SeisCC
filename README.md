# [Cross-correlation](https://en.wikipedia.org/wiki/Cross-correlation)

## cross-correlation methods

- time domain correlation
    - [cc_time](src/cc_time.c): only calculate cross-correlation coefficient at the specific time window

- frequncy domain correlation
    - [cc_freq](src/cc_freq.c): do cross-correlation using [FFTW](http://www.fftw.org/)
    - [correlatec](src/correlatec.c): call SAC's `crscor` in SAC's [libsac](http://ds.iris.edu/files/sac-manual/manual/saclib.html) ([Chinese notes](https://seisman.github.io/SAC_Docs_zh/libs/libsac/#crscor)) to do cross-correlation
    - [sac_wfcc](sac_wfcc/): this is [Lupei Zhu's code](http://geophysics.eas.gatech.edu/people/zpeng/Teaching/SAC_Tutorial/#part3_1)
    - PCC2: this is [Martin Schimmel's code](http://diapiro.ictja.csic.es/gt/mschi/SCIENCE/pcc2_method.py)
    - ObsPy: [corss_correlation](https://docs.obspy.org/packages/autogen/obspy.signal.cross_correlation.html#module-obspy.signal.cross_correlation)
    - SAC: [correlate](examples/SAC-correlate.sh)


### `cc_time`

```
Do correlation in time domain.

Usage:
  cc_time -Ttmark/ts/te [-h] sacifle1 sacfile2

Options:
  -T: tmark/begin time (sec)/time window (sec)
  -h: show usage
```

### `cc_freq`

```
Do correlation in frequency domain.

Usage:
  cc_freq -Ttmark/ts/te -Occf -Wtaper -Aautocorr
          [-h] sacfile1 sacfile2

Options:
  -T: tmark/begin time (sec)/time window (sec)
  -O: cross-correlation function file
  -W: taper (0: NO; 1: hanning; 2: cos)
  -A: output auto-correlation
  -h: show usage
```


### `correlatec`
```
Do correlation in frequency domain calling crscor.

Usage:
  correlatec -Ttmark/ts/te -Occf -Wtaper -Aautocorr
             [-h] sacfile1 sacfile2

Options:
  -T: tmark/begin time (sec)/time window (sec)
  -O: cross-correlation function file
  -W: taper (0: NO; 1: hanning; 2: cos)
  -A: output auto-correlation
  -h: show usage
```


## notes

When we do correlation, some pre-processing may have been done, e.g., `rmean`, `rtrend`, `taper`.

```bash
$ SAC
$ r seis.sac
$ rmean; rtrend; taper;
```
Some codes (e.g., `cc-time`) may just use the RAW data, while some codes (e.g., `cc_freq` and `correlatec`) may do some processing in the running, e.g. `taper`. Therefore, we have to be very careful to explain the difference between different codes if the data preprocessing is a little different. Maybe we should add those pre-processing in the code, so that we can choose whether to use it.

Be careful about reference time between different codes!!! Usually, we use the tmark (-5 -> b; -4 -> e; -3 -> 0; -2 -> a; 0-9 -> T0-9) as the reference time. But some code may add additional shift for some purpose, e.g., `sac_wfcc` adds tmark to the time shift.

