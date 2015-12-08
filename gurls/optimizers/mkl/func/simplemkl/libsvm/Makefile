#
# choose between intel compiler and gcc 
# 
CC = gcc-4.2
CFLAGS=-O3 -ffast-math -fomit-frame-pointer -fPIC #-mtune=k8 -march=k8

#CC=/agbs/share/sw/icc/bin/icc
#CFLAGS= -xP -fast -fPIC

MATLABDIR=/agbs/share/sw/matlab

MATLAB_INCLUDE=-I$(MATLABDIR)/extern/include
MATLAB_LIB=-L$(MATLABDIR)/bin/glnx86 -lmex -lm
#MATLAB_LIB=-L$(MATLABDIR)/bin/glnxa64 -lmex -lm

TARGETS=libsvm_classifier

INCLUDES=$(MATLAB_INCLUDE) -I.

# 32 bit version
SUFFIX=mexglx

# 64 bit version
#SUFFIX=mexa64


LDIRS=$(MATLAB_LIB) 

all: $(TARGETS:%=%.$(SUFFIX))

svm.o: svm.cpp svm.h
	$(CC) $(CFLAGS)  -c -o svm.o svm.cpp

%.o:	%.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) $(VERBOSEFLAG) -c -o $@ $*.cpp

%.$(SUFFIX):	svm.o mexarg.o svm.h mexarg.h %.o
	$(CC) $*.o $(LDIRS) -shared -o $*.$(SUFFIX) svm.o mexarg.o

### Targets

matlab_test:	libsvm_classifier.$(SUFFIX)
	matlab -nojvm -nodisplay -r "test_example; exit"; \
	matlab -nojvm -nodisplay -r "test_example_mc; exit";

clean:	
	rm -f *.mexglx *.mexa64 *.o *~

