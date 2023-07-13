DEF = -DDBL 
CFLAGS = -g $(DEF)

.c.o:
	gcc -c $(CFLAGS) $*.c

fir2caus:   fir2caus.o realft.o four1.o
	gcc $(CFLAGS) -o fir2caus fir2caus.o realft.o four1.o -lm

interpolate:   interpolate.o ipolstuff.o
	gcc $(CFLAGS) -o interpolate interpolate.o ipolstuff.o -lm
	
decimate:   decimate.o decim.o
	gcc $(CFLAGS) -o decimate decimate.o decim.o -lm

all: fir2caus interpolate decimate
