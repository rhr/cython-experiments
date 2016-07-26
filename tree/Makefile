CC	= gcc
CFLAGS	= -O3 -fPIC
FC	= gfortran
FFLAGS	= -O3 -fPIC

FOBJS	= expokit.o mataid.o clock.o
COBJS	= dexpm_c.o
CYSRC   = cyexpokit.c odeiv.c tree.c mcmc.c
CYMOD	= cyexpokit.so odeiv.so tree.so mcmc.so

%.c: %.pyx
	cython $<

%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $<

%.o : %.f90
	$(CC) $(CFLAGS) -c -o $@ $<

all: $(CYSRC) $(FOBJS) $(COBJS)
	python ./setup.py build_ext --inplace

clean: 
	rm $(FOBJS) $(COBJS) $(CYSRC) $(CYMOD)
