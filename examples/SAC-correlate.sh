#!/bin/bash
#
# calculate cross-correlation using correlate in SAC:
#   https://seisman.github.io/SAC_Docs_zh/commands/correlate/
#
# NOTES:
#   the final lag time (dT) is with respect to begin time, i.e.,
#      dT = dt + (b(f2) - b(f1))
#
# run the bash script to get arrival time
#  ./SAC-correlate.sh | awk '{print $2}'
#

# auto-correlation
sac << EOF
fg seis
w f1.sac
w f2.sac
r f1.sac f2.sac
echo on processed
cor norm
(max &2,depmax (abs &2,depmin))
q
EOF

rm f1.sac f2.sac


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


sac << EOF
r seis.sac seis-shift0.2345.sac
echo on processed
cor norm
(max &2,depmax (abs &2,depmin))
r seis-rrt.sac seis-shift0.2345-rrt.sac
cor norm
(max &2,depmax (abs &2,depmin))
q
EOF

rm seis*.sac

