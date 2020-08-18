
CC   = gcc
SUBS = sacio.o fft.o Complex.o

all: sac_wfcc src_ss clean

src_ss: src_ss.o $(SUBS)
	${CC} -o $@ $^ -lm

sac_wfcc: sac_wfcc.o $(SUBS)
	${CC} -o ../bin/$@ $^ -lm

clean:
	rm -f *.o

