
INCLUDESDIR = ./include

SRCDIR = ./src

OBJ = ./obj
NVCC = /usr/local/cuda/bin/nvcc

COMPILERGPU = $(NVCC)
COMPILERCPU = $(CC)


ifeq ($(IFT_GPU), 1)
    SRC = $(SRC_GPU)
    INCLUDE = $(INCLUDESDIR)  -I /usr/local/cuda/include
    FILES =  $(OBJ)/check_cuda_device.cu.o $(OBJ)/kernels.cu.o $(OBJ)/sptree.cu.o $(OBJ)/t_sne_gpu.cu.o
    MYFLAGS = -Wno-deprecated-gpu-targets -std=c++11
else
    SRC = $(SRC_CPU)
    INCLUDE = $(INCLUDESDIR)
    FILES = $(OBJ)/sptree.o $(OBJ)/tsne.o
    MYFLAGS = -fPIC -O2 -Wall
endif

all: $(FILES)

$(OBJ)/sptree.o: $(SRCDIR)/sptree.cpp
	$(COMPILERCPU) $(MYFLAGS) -c $< -I$(INCLUDE) -o $@

$(OBJ)/tsne.o: $(SRCDIR)/tsne.cpp
	$(COMPILERCPU) $(MYFLAGS) -c $< -I$(INCLUDE) -o $@

$(OBJ)/check_cuda_device.cu.o: $(SRCDIR)/check_cuda_device.cu
	$(COMPILERGPU) $(MYFLAGS) -c $< -I$(INCLUDE) -o $(OBJ)/check_cuda_device.o

$(OBJ)/kernels.cu.o: $(SRCDIR)/kernels.cu
	$(COMPILERGPU) $(MYFLAGS) -c $< -I$(INCLUDE) -o $(OBJ)/kernels.o

$(OBJ)/sptree.cu.o: $(SRCDIR)/sptree.cu
	$(COMPILERGPU) $(MYFLAGS) -c $< -I$(INCLUDE) -o $(OBJ)/sptree.o

$(OBJ)/t_sne_gpu.cu.o: $(SRCDIR)/t_sne_gpu.cu
	$(COMPILERGPU) $(MYFLAGS) -c $< -I$(INCLUDE) -o $(OBJ)/t_sne_gpu.o

clean:
	rm -f $(OBJ)/*.o;
