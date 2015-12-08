#
# choose between intel compiler and gcc 
# 

#CC = gcc-4.1
#CFLAGS=-O3 -ffast-math -fomit-frame-pointer -fPIC -mtune=k8 -march=k8

CC=/agbs/share/sw/icc/bin/icc
CFLAGS= -fast -Wall #-Werror

# your matlab dir here:
MATLABDIR=/agbs/share/sw/matlab 

# 32 bit version
SUFFIX=mexglx
LDIRS=-L$(MATLABDIR)/bin/glnx86 

# 64 bit version
#SUFFIX=mexa64
#LDIRS=-L$(MATLABDIR)/bin/glnxa64

TARGETS=Kbeta subproblem_product_objective
LIBS=-lm -lmex
INCLUDES=-I$(MATLABDIR)/extern/include -I.

all: $(TARGETS:%=%.$(SUFFIX)) 

%.o:	%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $*.cpp

%.o:	%.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $*.c

%.$(SUFFIX): %.o
	$(CC) $*.o $(LDIRS) $(LIBS) -shared -o $*.$(SUFFIX) 

clean:	
	rm -f *.mexglx *.mexa64 *.o *~ *.so

Kbeta_test:
	matlab -nojvm -nodisplay -r "addpath tests; Kbeta_test;exit"	

subprob_test:
	matlab -nojvm -nodisplay -r "addpath tests; subproblem_product_objective_test; exit"

libsvm_test: 
	cd libsvm; matlab -nojvm -nodisplay -r "test_example; exit";

ipopt_test:
	matlab -nojvm -nodisplay -r "addpath tests; ipopt_test;exit"

test: libsvm_test Kbeta_test subprob_test ipopt_test
