autocorr: autocorr.c autocorr.h
	cc autocorr.c -o autocorr -I . -lm -g3 -Wall

atuocorr.o: autocorr.c autocorr.h
	cc -c autocorr.c -I .

libatuocorr.a: autocorr.o
	ar -cr libautocorr.a autocorr.o

