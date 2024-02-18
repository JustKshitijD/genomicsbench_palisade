CC      = gcc
#CC		= icpc
ARCH_FLAGS=	-msse4.1

INC_PALISADE = -I/usr/local/include/palisade -I$(IN)/usr/local/include/palisade/cereal -I$(IN)/usr/local/include/palisade/pke -I$(IN)/usr/local/include/palisade/core
OPT_PALISADE = -DPALISADE_VERSION=1.11.7 -fopenmp
LIBS_PALISADE = -L/use/local/lib -lPALISADEbinfhe.so.1.11.7 -lPALISADEpke.so.1.11.7 -lPALISADEcore.so.1.11.7

ifeq ($(arch),sse41)
	ARCH_FLAGS=-msse4.1
else ifeq ($(arch),avx2)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-march=core-avx2 #-xCORE-AVX2
	else	
		ARCH_FLAGS=-mavx2
	endif
else ifeq ($(arch),avx512)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-xCORE-AVX512
	else	
		ARCH_FLAGS=-mavx512bw
	endif
else ifeq ($(arch),native)
		ARCH_FLAGS=-march=native
else ifneq ($(arch),)
	## To provide a different architecture flag like -march=core-avx2.
		ARCH_FLAGS=$(arch)
endif
CFLAGS= -g -O2 -std=c++11 $(ARCH_FLAGS)

#VTUNE_HOME= /opt/intel/vtune_profiler
GKLPATH= ../../tools/GKL
INC=
SHARED_LIBRARIES = -fopenmp -L$(GKLPATH)/build/native -lgkl_pairhmm_c -lstdc++ -ldl

ifneq ($(VTUNE_HOME),)
CFLAGS+= -DVTUNE_ANALYSIS=1
INC+= -I${VTUNE_HOME}/include
SHARED_LIBRARIES+= -L${VTUNE_HOME}/lib64 -littnotify
endif

all: phmm

phmm: PairHMMUnitTest.o /home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/palisade_header.o
	$(CC) $(INC_PALISADE) $(OPT_PALISADE) -o $@ $^ $(CFLAGS) ${INC} $(SHARED_LIBRARIES) /usr/local/lib/libPALISADEbinfhe.so.1.11.7 /usr/local/lib/libPALISADEpke.so.1.11.7 /usr/local/lib/libPALISADEcore.so.1.11.7

PairHMMUnitTest.o: PairHMMUnitTest.cpp PairHMMUnitTest.h
	$(CC) $(INC_PALISADE) $(OPT_PALISADE)  -c $(CFLAGS) ${INC} -fopenmp $< 

.PHONY: clean

clean:
	rm -f phmm *.o
