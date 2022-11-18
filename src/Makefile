
CC     = gcc
CFLAGS = -Wall

LIB    = -lfftw3 -lm

BIN    = ../bin


all: cc_time cc_freq correlatec clean

cc_time: cc_time.o sacio.o
	$(CC) ${CFLAGS} -o $(BIN)/$@ $^ -lm

cc_freq: cc_freq.o sacio.o taper.o
	$(CC) ${CFLAGS} -o $(BIN)/$@ $^ ${LIB}

correlatec: correlatec.o sacio.o taper.o
	$(CC) ${CFLAGS} -o $(BIN)/$@ $^ -lm -L/home/coreman/opt/sac-102/lib -lsac -lsacio

clean:
	rm *.o
