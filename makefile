##############################################################
# CUDA Sample Makefile                                       #
#                                                            #
# Author      : Jieqing Gan                                 #
# Version     : 3.5                                          #
# Date        : 10/01/2014                                   #
# Discription : generic Makefile for making CUDA programs    #
##############################################################

BIN               := ../GDEM


# Compilers
CUDA_PATH ?= /usr/local/cuda/8.0.61
MPI_PATH     ?= /usr/local/openmpi/1.10.7-mlx
NVCC       := $(CUDA_PATH)/bin/nvcc 
#CXX        := gcc -Xlinker --allow-multiple-definition
CXX        := $(MPI_PATH)/bin/mpic++ -Xlinker --allow-multiple-definition

MPICC       :=  $(MPI_PATH)/bin/mpic++ -Xlinker --allow-multiple-definition
EXEC   ?=

# Paths
INCD =  -I. -I$(CUDA_PATH)/include \
            -I$(MPI_PATH)/include
LIBS = -lm -L"$(CUDA_PATH)/lib64" -lcudart -g3

#-I/usr/lib64/openmpi/include

# internal flags
#NVCCFLAGS  :=
CCFLAGS     := $(INCD) -g3
NVCCLDFLAGS :=
LDFLAGS     :=

###########################################################

# CUDA code generation flags
GENCODE_FLAGS := -arch=compute_60 -code=sm_60,compute_60
#
# Program-specific
CFILES       := AllocateArrays.cpp \
                     Boundary.cpp  \
                     Feed.cpp  \
                     FreeArrays.cpp  \
                     GPUSet.cpp  \
                     Initialization.cpp  \
                     Materials.cpp  \
                     MemcpyHostDevice.cpp \
                     Movement.cpp  \
                                          mpiFunctions.cpp\
                     Particle.cpp  \
                     ReadData.cpp  \
                     WriteData.cpp  \
                     dempacking.cpp

CUFILES       := dempacking.cu

OFILES=$(CFILES:.cpp=.o)
CUOFILES=$(CUFILES:.cu=.cu_o)
USERCUOFILES= treatmesh.cu_o
##

OBJ  = $(OFILES) $(CUOFILES) $(USERCUOFILES)

# Build Rules
all: $(BIN) clean

$(BIN): $(OBJ)
	$(CXX) $(OBJ) -o  $(BIN) $(INCD) $(LIBS)

%.o : %.cpp 
	$(CXX) -c  $(CCFLAGS) -o $@ $<

%.cu_o : %.cu 
	$(NVCC) $(GENCODE_FLAGS) -c $(INCD) -o $@ $<

clean:
	rm -f $(OFILES) $(CUOFILES) *.o 
