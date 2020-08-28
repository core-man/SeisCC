#!/bin/bash
#
# compare time-domain and frequency-domain correlation codes
#

# prepare SAC data
# seis.sac, its shifted seis-shift[xxx].sac, and
# those two after do rmean, rtrend, taper
sac << EOF
fg seis
w seis.sac
r ./seis.sac
ch b 0
wh
cut 0 7
r seis.sac
w over
cut off
r seis.sac
ch t1 3
wh

cuterr fillz
cut 0.2345 7.2345
r seis.sac
w seis-shift0.2345.sac
cut off
r seis-shift0.2345.sac
ch b 0
w over
q
EOF

sac << EOF
r seis.sac seis-shift0.2345.sac
rmean
rtrend
taper
w seis-rrt.sac seis-shift0.2345-rrt.sac
q
EOF


# SAC's frequency-domain correlate
sac << EOF
r seis.sac seis-shift0.2345.sac
cor norm
w junk.sac SAC-ccf.sac
r seis-rrt.sac seis-shift0.2345-rrt.sac
cor norm
w junk.sac SAC-ccf-rrt.sac
q
EOF



# do time-domain and frequency cross-corrleation
## cross-correlation at zero lag time
### data without rmean;rtrend;taper
echo "----------------------------------"
echo "cross-correlation at zero lag time";
echo "----------------------------------"
echo "data without rmean;rtrend;taper:"
echo "1. cc_time: time-domain cross-correlation at zero lag time (f1 f2 ccf)"
../bin/cc_time -T1/-3/4 seis.sac seis-shift0.2345.sac
echo
echo "2. cc_freq: frequency-domain cross-correlation at zero lag time (f1 f2 ccf)"
../bin/cc_freq -T1/-3/4 -A0 seis.sac seis-shift0.2345.sac
echo
echo "3. correlatec: frequency-domain cross-correlation at zero lag time (f1 f2 ccf)"
../bin/correlatec -T1/-3/4 -A0 seis.sac seis-shift0.2345.sac
echo
echo "4. SAC-correlate: frequency-domain cross-correlation at zero lag time (cross-correlation-function ccf)"
/home/tomoboy/src.import/SACTools/bin/sacamp -T0 SAC-ccf.sac

### data with rmean;rtrend;taper
echo
echo "data with rmean;rtrend;taper:"
echo "1. cc_time: time-domain cross-correlation at zero lag time (f1 f2 ccf)"
../bin/cc_time -T1/-3/4 seis-rrt.sac seis-shift0.2345-rrt.sac
echo
echo "2. cc_freq: frequency-domain cross-correlation at zero lag time (f1 f2 ccf)"
../bin/cc_freq -T1/-3/4 -A0 seis-rrt.sac seis-shift0.2345-rrt.sac
echo
echo "3. correlatec: frequency-domain cross-correlation at zero lag time (f1 f2 ccf)"
../bin/correlatec -T1/-3/4 -A0 seis-rrt.sac seis-shift0.2345-rrt.sac
echo
echo "4. SAC-correlate: frequency-domain cross-correlation at zero lag time (cross-correlation-function ccf)"
/home/tomoboy/src.import/SACTools/bin/sacamp -T0 SAC-ccf-rrt.sac


## cross-correlation for all the lag times
### data without rmean;rtrend;taper
echo
echo
echo "-----------------------------------"
echo "cross-correlation all the lag times";
echo "-----------------------------------"
echo "data without rmean;rtrend;taper:"
echo "1. cc_freq: frequency-domain cross-correlation (f1 f2 ccf dt)"
../bin/cc_freq -T1/-3/4 seis.sac seis-shift0.2345.sac
echo
echo "2. correlatec: frequency-domain cross-correlation (f1 f2 ccf dt)"
../bin/correlatec -T1/-3/4 seis.sac seis-shift0.2345.sac
echo
echo "3. SAC-correlate: frequency-domain cross-correlation at -0.23 lag time (cross-correlation-function ccf)"
/home/tomoboy/src.import/SACTools/bin/sacamp -T-0.23 SAC-ccf.sac
echo
echo "4. Lupei Zhu's sac_wfcc"
echo "   f1  tmark     ccf"
echo "   f2  tmark+dt  ccf"
saclst t1 f seis.sac seis-shift0.2345.sac | ../bin/sac_wfcc -D-3/4/4

### data with rmean;rtrend;taper
echo
echo "data with rmean;rtrend;taper"
echo "1. cc_freq: frequency-domain cross-correlation (f1 f2 ccf dt)"
../bin/cc_freq -T1/-3/4 seis-rrt.sac seis-shift0.2345-rrt.sac
echo
echo "2. correlatec: frequency-domain cross-correlation (f1 f2 ccf dt)"
../bin/correlatec -T1/-3/4 seis-rrt.sac seis-shift0.2345-rrt.sac
echo
echo "3. SAC-correlate: frequency-domain cross-correlation at -0.23 lag time (cross-correlation-function ccf)"
/home/tomoboy/src.import/SACTools/bin/sacamp -T-0.23 SAC-ccf-rrt.sac
echo
echo "4. Lupei Zhu's sac_wfcc"
echo "   f1  tmark     ccf"
echo "   f2  tmark+dt  ccf"
saclst t1 f seis-rrt.sac seis-shift0.2345-rrt.sac | ../bin/sac_wfcc -D-3/4/4


rm *.sac

