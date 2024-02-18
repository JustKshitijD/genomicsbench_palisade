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

INC_PALISADE = -I/usr/local/include/palisade -I$(IN)/usr/local/include/palisade/cereal -I$(IN)/usr/local/include/palisade/pke -I$(IN)/usr/local/include/palisade/core
OPT_PALISADE = -DPALISADE_VERSION=1.11.7 -fopenmp
LIBS_PALISADE = -L/use/local/lib -lPALISADEbinfhe.so.1.11.7 -lPALISADEpke.so.1.11.7 -lPALISADEcore.so.1.11.7

.PHONY:all
	
all: create_packed_CT

create_packed_CT: palisade_header.o
	g++ -g $(INC_PALISADE) $(OPT_PALISADE)  create_packed_CT.cpp -o $@ $< /usr/local/lib/libPALISADEbinfhe.so.1.11.7 /usr/local/lib/libPALISADEpke.so.1.11.7 /usr/local/lib/libPALISADEcore.so.1.11.7

palisade_header.o: palisade_header.cpp palisade_header.h 
	g++ -g -c $(INC_PALISADE) $(OPT_PALISADE)  $< -o $@

# palisade_header.o: palisade_header.cpp
# 	g++ -c -O2 $(INC_PALISADE) $(OPT_PALISADE) palisade_header.cpp -o  $@ 
