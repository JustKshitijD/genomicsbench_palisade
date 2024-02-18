CXX=g++
#CXX=icpc

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

CXXFLAGS=-O3 -fopenmp $(ARCH_FLAGS) -g
#VTUNE_HOME= /opt/intel/vtune_profiler
INCLUDES=-I../../tools/spoa/include -I.
LIBS=-L../../tools/spoa/build/lib/ -lspoa -ldl -fopenmp

ifneq ($(VTUNE_HOME),)
CXXFLAGS+= -DVTUNE_ANALYSIS=1
INCLUDES+=	-I${VTUNE_HOME}/include
LIBS+=-L${VTUNE_HOME}/lib64 -littnotify
endif

all: msa_spoa_omp.cpp palisade_header.o
	$(CXX) $(INC_PALISADE) $(OPT_PALISADE) $(CXXFLAGS) msa_spoa_omp.cpp ${INCLUDES} ${LIBS} -o poa palisade_header.o /usr/local/lib/libPALISADEbinfhe.so.1.11.7 /usr/local/lib/libPALISADEpke.so.1.11.7 /usr/local/lib/libPALISADEcore.so.1.11.7

.PHONY: clean

clean:
	rm -f poa

palisade_header.o: palisade_header.cpp
	$(info compiling palisade_header.o)
	$(CXX) $(INC_PALISADE) $(OPT_PALISADE) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) palisade_header.cpp -o palisade_header.o
	$(info done compiling)

