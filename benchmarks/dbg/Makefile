CXX     = g++
#CXX		 = icpc
#VTUNE_HOME= /opt/intel/vtune_profiler
LIBS = -L../../tools/htslib -lcurl -lcrypto -lz -lm -lbz2 -llzma -lpthread -ldl -fopenmp
INC = -I../../tools/htslib -I.
CXXFLAGS   = -g -Wall -O2 -fopenmp

INC_PALISADE = -I/usr/local/include/palisade -I$(IN)/usr/local/include/palisade/cereal -I$(IN)/usr/local/include/palisade/pke -I$(IN)/usr/local/include/palisade/core
OPT_PALISADE = -DPALISADE_VERSION=1.11.7 -fopenmp
LIBS_PALISADE = -L/use/local/lib -lPALISADEbinfhe.so.1.11.7 -lPALISADEpke.so.1.11.7 -lPALISADEcore.so.1.11.7

ifneq ($(VTUNE_HOME),)
CXXFLAGS+= -DVTUNE_ANALYSIS=1
INC+=	-I${VTUNE_HOME}/include
LIBS+=-L${VTUNE_HOME}/lib64 -littnotify
endif

all: debruijn.cpp common.cpp common.h palisade_header.o
	$(CXX) debruijn.cpp common.cpp ../../tools/htslib/libhts.a $(CXXFLAGS) $(INC) $(LIBS) -o dbg palisade_header.o /usr/local/lib/libPALISADEbinfhe.so.1.11.7 /usr/local/lib/libPALISADEpke.so.1.11.7 /usr/local/lib/libPALISADEcore.so.1.11.7

clean: 
	rm -f dbg

palisade_header.o: palisade_header.cpp
	$(info compiling palisade_header.o)
	$(CXX) $(INC_PALISADE) $(OPT_PALISADE) -c $(CXXFLAGS) $(CPPFLAGS) $(INC) palisade_header.cpp -o palisade_header.o
	$(info done compiling)

