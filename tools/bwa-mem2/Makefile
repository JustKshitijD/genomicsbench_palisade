##/*************************************************************************************
##                           The MIT License
##
##   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
##   Copyright (C) 2019  Intel Corporation, Heng Li.
##
##   Permission is hereby granted, free of charge, to any person obtaining
##   a copy of this software and associated documentation files (the
##   "Software"), to deal in the Software without restriction, including
##   without limitation the rights to use, copy, modify, merge, publish,
##   distribute, sublicense, and/or sell copies of the Software, and to
##   permit persons to whom the Software is furnished to do so, subject to
##   the following conditions:
##
##   The above copyright notice and this permission notice shall be
##   included in all copies or substantial portions of the Software.
##
##   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
##   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
##   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
##   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
##   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
##   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
##   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
##   SOFTWARE.
##
##Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
##                                Heng Li <hli@jimmy.harvard.edu> 
##*****************************************************************************************/

ifneq ($(portable),)
	STATIC_GCC=-static-libgcc -static-libstdc++
endif

INC_PALISADE = -I/usr/local/include/palisade -I$(IN)/usr/local/include/palisade/cereal -I$(IN)/usr/local/include/palisade/pke -I$(IN)/usr/local/include/palisade/core
OPT_PALISADE = -DPALISADE_VERSION=1.11.7 -fopenmp
LIBS_PALISADE = -L/usr/local/lib -lPALISADEbinfhe.so.1.11.7 -lPALISADEpke.so.1.11.7 -lPALISADEcore.so.1.11.7

EXE=		bwa-mem2
#CXX=		icpc
ifeq ($(CXX), icpc)
	CC= icc
else ifeq ($(CXX), g++)
	CC=gcc
endif		
ARCH_FLAGS=	-msse4.1
MEM_FLAGS=	-DSAIS=1
CPPFLAGS+=	-DENABLE_PREFETCH -DV17=1 $(MEM_FLAGS) 
INCLUDES=   -Isrc -Iext/safestringlib/include
LIBS=		-lpthread -lm -lz -L. -lbwa -Lext/safestringlib -lsafestring $(STATIC_GCC) -fopenmp
OBJS=		/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/palisade_header.o src/fastmap.o src/bwtindex.o src/utils.o src/memcpy_bwamem.o src/kthread.o \
			src/kstring.o src/ksw.o src/bntseq.o src/bwamem.o src/profiling.o src/bandedSWA.o \
			src/FMI_search.o src/read_index_ele.o src/bwamem_pair.o src/kswv.o src/bwa.o \
			src/bwamem_extra.o src/kopen.o
BWA_LIB=    libbwa.a						# include "/usr/local/include/palisade/cereal"
SAFE_STR_LIB=    ext/safestringlib/libsafestring.a

ifeq ($(arch),sse41)
	ARCH_FLAGS=-msse4.1
else ifeq ($(arch),sse42)
	ARCH_FLAGS=-msse4.2
else ifeq ($(arch),avx)
	ifeq ($(CXX), icpc)
		ARCH_FLAGS=-mavx ##-xAVX
	else	
		ARCH_FLAGS=-mavx
	endif
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
# To provide a different architecture flag like -march=core-avx2.
	ARCH_FLAGS=$(arch)
else

endif

CXXFLAGS+=	-g -O3 -fpermissive $(ARCH_FLAGS) #-Wall ##-xSSE2

.PHONY:all clean depend multi
.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) $(INC_PALISADE) $(OPT_PALISADE) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@ 
# /usr/local/lib/libPALISADEbinfhe.so.1.11.7 /usr/local/lib/libPALISADEpke.so.1.11.7 /usr/local/lib/libPALISADEcore.so.1.11.7

# The header files to the right means that if any of these header files change, the .o file will have to be compiled again
# The -c and -o $@ means that object file is going to be created (due to -c) and name of that will be $@(the left of :)
# $< is the first element to right of : , which is .cpp file in this case

all:$(EXE)

multi:
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse41    EXE=bwa-mem2.sse41    CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=sse42    EXE=bwa-mem2.sse42    CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx    EXE=bwa-mem2.avx    CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx2   EXE=bwa-mem2.avx2     CXX=$(CXX) all
	rm -f src/*.o $(BWA_LIB); cd ext/safestringlib/ && $(MAKE) clean;
	$(MAKE) arch=avx512 EXE=bwa-mem2.avx512bw CXX=$(CXX) all
	$(CXX) -Wall -O3 src/runsimd.cpp -Iext/safestringlib/include -Lext/safestringlib/ -lsafestring $(STATIC_GCC) -o bwa-mem2

$(EXE):$(BWA_LIB) $(SAFE_STR_LIB) src/main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) src/main.o $(BWA_LIB) $(LIBS) -o $@ /usr/local/lib/libPALISADEbinfhe.so.1.11.7 /usr/local/lib/libPALISADEpke.so.1.11.7 /usr/local/lib/libPALISADEcore.so.1.11.7

# ar -rcs is the most likely command you would use when using a Makefile to compile a library. 
# r means that if the library already exists, replace the old files within the library with your new files. c means create the library if it did not exist. 
# s can be seen to mean 'sort' (create a sorted index of) the library, so that it will be indexed and faster to access the functions in the library. 
# $(OBJS) is all the 18 object files
# BWA_LIB=    libbwa.a
$(BWA_LIB):$(OBJS)
	ar rcs $(BWA_LIB) $(OBJS) 

$(SAFE_STR_LIB):
	cd ext/safestringlib/ && $(MAKE) clean && $(MAKE) CC=$(CC) directories libsafestring.a

clean:
	rm -fr src/*.o $(BWA_LIB) $(EXE) bwa-mem2.sse41 bwa-mem2.avx2 bwa-mem2.avx512bw
	cd ext/safestringlib/ && $(MAKE) clean

depend:
	(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CXXFLAGS) $(CPPFLAGS) -I. -- src/*.cpp)

# DO NOT DELETE

src/FMI_search.o: src/FMI_search.h src/bntseq.h src/read_index_ele.h 
src/FMI_search.o: src/utils.h src/macro.h src/bwa.h src/bwt.h src/sais.h
src/bandedSWA.o: src/bandedSWA.h src/macro.h
src/bntseq.o: src/bntseq.h src/utils.h src/macro.h src/kseq.h src/khash.h src/kstring.h
src/bwa.o: src/bntseq.h src/bwa.h src/bwt.h src/macro.h src/ksw.h src/utils.h
src/bwa.o: src/kstring.h src/kvec.h src/kseq.h
src/bwamem.o: src/bwamem.h src/bwt.h src/bntseq.h src/bwa.h src/macro.h
src/bwamem.o: src/kthread.h src/bandedSWA.h src/kstring.h src/ksw.h
src/bwamem.o: src/kvec.h src/ksort.h src/utils.h src/profiling.h
src/bwamem.o: src/FMI_search.h src/read_index_ele.h src/kbtree.h
src/bwamem_extra.o: src/bwa.h src/bntseq.h src/bwt.h src/macro.h src/bwamem.h
src/bwamem_extra.o: src/kthread.h src/bandedSWA.h src/kstring.h src/ksw.h
src/bwamem_extra.o: src/kvec.h src/ksort.h src/utils.h src/profiling.h
src/bwamem_extra.o: src/FMI_search.h src/read_index_ele.h
src/bwamem_pair.o: src/kstring.h src/bwamem.h src/bwt.h src/bntseq.h
src/bwamem_pair.o: src/bwa.h src/macro.h src/kthread.h src/bandedSWA.h
src/bwamem_pair.o: src/ksw.h src/kvec.h src/ksort.h src/utils.h
src/bwamem_pair.o: src/profiling.h src/FMI_search.h src/read_index_ele.h
src/bwamem_pair.o: src/kswv.h
src/bwtindex.o: src/bntseq.h src/bwa.h src/bwt.h src/macro.h src/utils.h
src/bwtindex.o: src/FMI_search.h src/read_index_ele.h
src/fastmap.o: src/fastmap.h src/bwa.h src/bntseq.h src/bwt.h src/macro.h
src/fastmap.o: src/bwamem.h src/kthread.h src/bandedSWA.h src/kstring.h
src/fastmap.o: src/ksw.h src/kvec.h src/ksort.h src/utils.h src/profiling.h
src/fastmap.o: src/FMI_search.h src/read_index_ele.h src/kseq.h src/kstring.h
src/kstring.o: src/kstring.h
src/ksw.o: src/ksw.h src/macro.h
src/kswv.o: src/kswv.h src/macro.h src/ksw.h src/bandedSWA.h
src/kthread.o: src/kthread.h src/macro.h src/bwamem.h src/bwt.h src/bntseq.h
src/kthread.o: src/bwa.h src/bandedSWA.h src/kstring.h src/ksw.h src/kvec.h
src/kthread.o: src/ksort.h src/utils.h src/profiling.h src/FMI_search.h
src/kthread.o: src/read_index_ele.h
src/main.o: src/main.h src/kstring.h src/utils.h src/macro.h src/bandedSWA.h
src/main.o: src/profiling.h
src/profiling.o: src/macro.h
src/read_index_ele.o: src/read_index_ele.h src/utils.h src/bntseq.h
src/read_index_ele.o: src/macro.h
src/utils.o: src/utils.h src/ksort.h src/kseq.h src/kstring.h
src/memcpy_bwamem.o: src/memcpy_bwamem.h
/home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/palisade_header.o: /home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/palisade_header.h
