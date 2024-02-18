# Run via "make CXX=g++ arch=avx2"

#CXX = icpc
CXX = g++

# path #
SRC_PATH = src
BUILD_PATH = build
BIN_PATH = $(BUILD_PATH)/bin

INC_PALISADE = -I/usr/local/include/palisade -I$(IN)/usr/local/include/palisade/cereal -I$(IN)/usr/local/include/palisade/pke -I$(IN)/usr/local/include/palisade/core
OPT_PALISADE = -DPALISADE_VERSION=1.11.7 -fopenmp
LIBS_PALISADE = -L/use/local/lib -lPALISADEbinfhe.so.1.11.7 -lPALISADEpke.so.1.11.7 -lPALISADEcore.so.1.11.7

# executable #
BIN_NAME = chain

# extensions #
SRC_EXT = cpp

# code lists #
# Find all source files in the source directory, sorted by
# most recently modified
SOURCES = $(shell find $(SRC_PATH) -name '*.$(SRC_EXT)' | sort -k 1nr | cut -f2-)
# SOURCES += /home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/benchmarks/chain/palisade_header.cpp
# SOURCES	+= /home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/palisade_header.cpp
# Set the object file names, with the source directory stripped
# from the path, and the build path prepended in its place

OBJECTS = $(SOURCES:$(SRC_PATH)/%.$(SRC_EXT)=$(BUILD_PATH)/%.o)
# OBJECTS += /home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/benchmarks/chain/palisade_header.o
# OBJECTS = build/main.o build/common.o build/host_data_io.o build/host_kernel.o /home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/palisade_header.o

# Set the dependency files that will be used to add header dependencies
DEPS = $(OBJECTS:.o=.d)


ARCH_FLAGS=	-msse4.1

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

# flags #
#COMPILE_FLAGS = -std=c++11 -Wall -Wextra -g -O3 -fopenmp -xAVX2 -axAVX2
#VTUNE_HOME= /opt/intel/vtune_profiler
COMPILE_FLAGS = -std=c++11 -Wall -Wextra -g -O3 -fopenmp $(ARCH_FLAGS)
INCLUDES = -I../../tools/minimap2 -I.
# Space-separated pkg-config libraries used by this project
LIBS = -L../../tools/minimap2 -lminimap2 -ldl


ifneq ($(VTUNE_HOME),)
COMPILE_FLAGS += -DVTUNE_ANALYSIS=1
INCLUDES += -I${VTUNE_HOME}/include
LIBS +=  -L${VTUNE_HOME}/lib64 -littnotify
endif


.PHONY: default_target
default_target: release

.PHONY: release
release: export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS)
release: dirs
	@$(MAKE) all

.PHONY: dirs
dirs:
	@echo "Creating directories"
	@mkdir -p $(dir $(OBJECTS))
	@mkdir -p $(BIN_PATH)

.PHONY: clean
clean:
	@echo "Deleting $(BIN_NAME) symlink"
	@$(RM) $(BIN_NAME)
	@echo "Deleting directories"
	@$(RM) -r $(BUILD_PATH)
	@$(RM) -r $(BIN_PATH)

# checks the executable and symlinks to the output
# $(BIN_NAME): chain; $(BIN_PATH)/$(BIN_NAME)=build/bin/chain
.PHONY: all
all: $(BIN_PATH)/$(BIN_NAME)
	@echo "Making symlink: $(BIN_NAME) -> $<"
	@$(RM) $(BIN_NAME)
	@ln -s $(BIN_PATH)/$(BIN_NAME) $(BIN_NAME)

# Creation of the executable
$(BIN_PATH)/$(BIN_NAME): $(OBJECTS)
	@echo "Linking: $@"
	#$(CXX) -O3 -fopenmp -xAVX2 -axAVX2 $(OBJECTS) -o $@
	$(CXX) -O3 -fopenmp $(ARCH_FLAGS) $(OBJECTS) $(INCLUDES) $(LIBS) -o $@ /usr/local/lib/libPALISADEbinfhe.so.1.11.7 /usr/local/lib/libPALISADEpke.so.1.11.7 /usr/local/lib/libPALISADEcore.so.1.11.7

# Add dependency files, if they exist
-include $(DEPS)

# Source file rules
# After the first compilation they will be joined with the rules from the
# dependency files to provide header dependencies

# OBJ_1 = $(BUILD_PATH)/%.o
# OBJ_1 += /home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/palisade_header.o
# SRC_1 = $(SRC_PATH)/%.$(SRC_EXT)	
# SRC_1 += /home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/palisade_header.cpp

$(BUILD_PATH)/%.o: $(SRC_PATH)/%.$(SRC_EXT)
	@echo "Compiling: $< -> $@"
	$(CXX) $(INC_PALISADE) $(OPT_PALISADE)  $(CXXFLAGS) $(INCLUDES) $(LIBS) -MP -MMD -c $< -o $@

# build/main.o: src/host_data_io.h src/host_data.h src/host_kernel.h
# build/common.o: src/common.h
# build/host_data_io.o: src/host_data_io.h src/host_data.h
# build/host_kernel.o: src/host_kernel.h src/common.h
# /home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/palisade_header.o: /home/kshitij/Desktop/NTU/Homomorphic/genomicsbench_palisade/palisade_header.h

# Flow-
# Compiling: src/common.cpp -> build/common.o
# Compiling: src/host_data_io.cpp -> build/host_data_io.o
# Compiling: src/host_kernel.cpp -> build/host_kernel.o
# Compiling: src/main.cpp -> build/main.o
# Linking: build/bin/chain
# Making symlink: chain -> build/bin/chain
