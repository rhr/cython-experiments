CC	= gcc
CFLAGS	= -O3 -fPIC
FC	= gfortran
FFLAGS	= -O3 -fPIC

LIBS	= -llapack -lblas
FOBJS	= matexp.o expokit.o mataid.o clock.o
CYMOD	= cyexpokit.so

all:	$(FOBJS) $(CYMOD)

clean: 
	rm *.o

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o : %.f90
	$(CC) $(CFLAGS) -c -o $@ $<

$(CYMOD): cyexpokit.pyx $(FOBJS) $(COBJS)
	rm -f cyexpokit.c
	python ./setup.py build_ext --inplace
