# [Cross-correlation](https://en.wikipedia.org/wiki/Cross-correlation)

## Cross-correlation methods

Time domain correlation

- [cc_time](src/cc_time.c): only calculate cross-correlation coefficient at the specific time window

Frequncy domain correlation

- [cc_freq](src/cc_freq.c): Do cross-correlation using [FFTW](http://www.fftw.org/)
- [correlatec](src/correlatec.c): Call the `crscor` function in SAC's [libsac](http://ds.iris.edu/files/sac-manual/manual/saclib.html) to do the cross-correlation. Please refer to [Chinese SAC Documentation](https://seisman.github.io/SAC_Docs_zh/libs/libsac/#crscor) for more details.
- [sac_wfcc](sac_wfcc/): This is [Lupei Zhu](http://www.eas.slu.edu/People/LZhu/home.html)'s code downloaded from [Introduction to SAC](http://geophysics.eas.gatech.edu/people/zpeng/Teaching/SAC_Tutorial/#part3_1) on 2020 Aug. 18.
- [PCC2](pcc/pcc2_method.py): This is [Martin Schimmel](http://diapiro.ictja.csic.es/gt/mschi/)'s code downloaded from [here](http://diapiro.ictja.csic.es/gt/mschi/SCIENCE/tseries.html#software) on 2021 Feb. 22.
- ObsPy: [corss_correlation](https://docs.obspy.org/packages/autogen/obspy.signal.cross_correlation.html#module-obspy.signal.cross_correlation)
- SAC: [correlate](examples/SAC-correlate.sh)

### `cc_time`

```
Do correlation in time domain.

Usage:
  cc_time -Ttmark/ts/te [-h] sacifle1 sacfile2

Options:
  -T: tmark/begin time (sec)/end time (sec)
  -h: show usage

Examples:
  cc_time -T1/-3/4 seis1.sac seis2.sac
```

### `cc_freq`

```
Do correlation in frequency domain.

Usage:
  cc_freq -Ttmark/ts/te [-Occf] [-Wtaper] [-Acczero]
          [-h] sacfile1 sacfile2

Options:
  -T: tmark/begin time (sec)/end time (sec)
  -O: cross-correlation function file
  -W: taper (0: NO (default); 1: hanning; 2: cos)
  -A: only output cross-correlation at zero lag time (0: YES; 1: NO (default))
  -h: show usage

Examples:
  cc_freq -T1/-3/4 seis1.sac seis2.sac
  cc_freq -T1/-3/4 -A0 seis1.sac seis2.sac
  cc_freq -T1/-3/4 -W1 seis1.sac seis2.sac
  cc_freq -T1/-3/4 -Occf12.sac seis1.sac seis2.sac
```


### `correlatec`

```
Do correlation in frequency domain calling crscor.

Usage:
  correlatec -Ttmark/ts/te [-Occf] [-Wtaper] [-Acczero]
             [-h] sacfile1 sacfile2

Options:
  -T: tmark/begin time (sec)/end time (sec)
  -O: cross-correlation function file
  -W: taper (0: NO (default); 1: hanning; 2: cos)
  -A: only output cross-correlation at zero lag time (0: YES; 1: NO (default))
  -h: show usage

Examples:
  correlatec -T1/-3/4 seis1.sac seis2.sac
  correlatec -T1/-3/4 -A0 seis1.sac seis2.sac
  correlatec -T1/-3/4 -W1 seis1.sac seis2.sac
  correlatec -T1/-3/4 -Occf12.sac seis1.sac seis2.sac
```


## Notes

When we do correlation, some pre-processing may have been done, e.g., `rmean`, `rtrend`, `taper`.

```bash
$ SAC
$ r seis.sac
$ rmean; rtrend; taper;
```
Some codes (e.g., `cc-time`) may just use the RAW data, while some codes (e.g., `cc_freq` and `correlatec`) may do some processing in the running, e.g. `taper`. Therefore, we have to be very careful to explain the difference between different codes if the data preprocessing is a little different. Maybe we should add those pre-processing in the code, so that we can choose whether to use it.

Be careful about reference time between different codes!!! Usually, we use the tmark (-5 -> b; -4 -> e; -3 -> 0; -2 -> a; 0-9 -> T0-9) as the reference time. But some code may add additional shift for some purpose, e.g., `sac_wfcc` adds tmark to the time shift.
