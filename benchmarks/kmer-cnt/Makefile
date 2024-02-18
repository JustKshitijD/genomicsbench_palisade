CXX     = g++
#CXX		 = icpc
#VTUNE_HOME= /opt/intel/vtune_profiler
LIBS= -lz -lm -lbz2 -llzma -ldl -fopenmp
CXXFLAGS   = -O3 -fopenmp
INC = -Ilibcuckoo

INC_PALISADE = -I/usr/local/include/palisade -I$(IN)/usr/local/include/palisade/cereal -I$(IN)/usr/local/include/palisade/pke -I$(IN)/usr/local/include/palisade/core
OPT_PALISADE = -DPALISADE_VERSION=1.11.7 -fopenmp

OBJS = palisade_header.o sequence_container.o sequence.o vertex_index.o 

FILE_LIB = libbwa.a

ifneq ($(VTUNE_HOME),)
CXXFLAGS+= -DVTUNE_ANALYSIS=1
LIBS += -L${VTUNE_HOME}/lib64 -littnotify
INC += -I${VTUNE_HOME}/include
endif

all: $(FILE_LIB) kmer_cnt.o
	$(CXX) $(CXXFLAGS) $(INC_PALISADE) $(OPT_PALISADE) kmer_cnt.o $(INC) $(FILE_LIB) $(LIBS) -o kmer-cnt /usr/local/lib/libPALISADEbinfhe.so.1.11.7 /usr/local/lib/libPALISADEpke.so.1.11.7 /usr/local/lib/libPALISADEcore.so.1.11.7

%.o:%.cpp
	$(info Hello)
	$(CXX)  $(CXXFLAGS) $(INC_PALISADE) $(OPT_PALISADE) $(INC) $(LIBS) -c $< -o $@

.PHONY: clean

clean: 
	rm -f kmer-cnt

$(FILE_LIB):$(OBJS)
	ar rcs $(FILE_LIB) $(OBJS)

palisade_header.o: palisade_header.h
sequence_container.o: palisade_header.h sequence_container.h sequence.h logger.h
sequence.o: palisade_header.h sequence.h
vertex_index.o: palisade_header.h vertex_index.h parallel.h memory_info.h logger.h config.h sequence_container.h kmer.h
kmer-cnt.o: palisade_header.h sequence_container.h sequence.h vertex_index.h memory_info.h logger.h utils.h config.h parallel.h
