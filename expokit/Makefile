CC	= gcc
CFLAGS	= -O3 -fPIC
FC	= gfortran
FFLAGS	= -O3 -fPIC

LIBS	= -llapack -lblas
FOBJS	= expokit.o mataid.o clock.o
COBJS	= dexpm_c.o
CYMOD	= cyexpokit.so tree.so

all:	$(CYMOD) $(FOBJS) $(COBJS)

clean: 
	rm *.o

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o : %.f90
	$(CC) $(CFLAGS) -c -o $@ $<

$(CYMOD): cyexpokit.pyx $(FOBJS) $(COBJS)
	rm -f cyexpokit.c
	python ./setup.py build_ext --inplace
