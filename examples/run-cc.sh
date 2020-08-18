#!/bin/bash
#
# compare time-domain and frequency-domain cc (SAC's correlate)
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



# time-domain and frequency cross-corrleation
echo
echo cc_time: time-domain cross-correlation at zero lag time
../bin/cc_time -T1/-3/4 seis.sac seis-shift0.2345.sac
echo cc_freq: frequency-domain cross-correlation at zero lag time
../bin/cc_freq -T1/-3/4 -W1 -A0 seis.sac seis-shift0.2345.sac
echo correlatec: frequency-domain cross-correlation at zero lag time
../bin/correlatec -T1/-3/4 -W1 -A0 seis.sac seis-shift0.2345.sac
echo SAC-correlate: frequency-domain cross-correlation at zero lag time
/home/tomoboy/src.import/SACTools/bin/sacamp -T0 SAC-ccf.sac

echo
echo cc_time: time-domain cross-correlation at zero lag time
../bin/cc_time -T1/-3/4 seis-rrt.sac seis-shift0.2345-rrt.sac
echo cc_freq: frequency-domain cross-correlation at zero lag time
../bin/cc_freq -T1/-3/4 -W1 -A0 seis-rrt.sac seis-shift0.2345-rrt.sac
echo correlatec: frequency-domain cross-correlation at zero lag time
../bin/correlatec -T1/-3/4 -W1 -A0 seis-rrt.sac seis-shift0.2345-rrt.sac
echo SAC-correlate: frequency-domain cross-correlation at zero lag time
/home/tomoboy/src.import/SACTools/bin/sacamp -T0 SAC-ccf-rrt.sac

echo ----------------------------------

echo
echo cc_freq: frequency-domain cross-correlation
../bin/cc_freq -T1/-3/4 -W1 -Occ_freq-ccf.sac  seis.sac seis-shift0.2345.sac
echo correlatec: frequency-domain cross-correlation
../bin/correlatec -T1/-3/4 -W1 -A1 -Ocorrelatec-ccf.sac seis.sac seis-shift0.2345.sac
echo SAC-correlate: frequency-domain cross-correlation at -0.23 lag time
/home/tomoboy/src.import/SACTools/bin/sacamp -T-0.23 SAC-ccf.sac
echo Lupei-sac_wfcc
saclst t1 f seis.sac seis-shift0.2345.sac | ../bin/sac_wfcc -D-3/4/4

echo
echo cc_freq: frequency-domain cross-correlation
../bin/cc_freq -T1/-3/4 -W1 -Occ_freq-ccf-rrt.sac seis-rrt.sac seis-shift0.2345-rrt.sac
echo correlatec: frequency-domain cross-correlation
../bin/correlatec -T1/-3/4 -W1 -A1 -Ocorrelatec-ccf.sac seis-rrt.sac seis-shift0.2345-rrt.sac
echo SAC-correlate: frequency-domain cross-correlation at -0.23 lag time
/home/tomoboy/src.import/SACTools/bin/sacamp -T-0.23 SAC-ccf-rrt.sac
echo Lupei-sac_wfcc
saclst t1 f seis-rrt.sac seis-shift0.2345-rrt.sac | ../bin/sac_wfcc -D-3/4/4


rm *.sac

