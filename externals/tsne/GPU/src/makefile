
CUDAPATH=/usr/local/cuda

# Have this point to an old enough gcc (for nvcc)
GCCPATH=/usr

NVCC=/usr/local/cuda/bin/nvcc
CCPATH=${GCCPATH}/bin

INCLUDES = -I/usr/local/cuda/include
SOURCES = nvmatrix.cu  nvmatrix_kernel.cu  tsne_p.cu 

#-O3 -lcuda -lcublas -lcudart 
FLAGS = -Wno-deprecated-gpu-targets -std=c++11

all: check_cuda_device.o kernels.o sptree.o t_sne_gpu.o

check_cuda_device.o: check_cuda_device.cu 
	$(NVCC)  $(INCLUDES) $(FLAGS) -c -o check_cuda_device.o check_cuda_device.cu

kernels.o: kernels.cu 
	$(NVCC)  $(INCLUDES) $(FLAGS) -c -o kernels.o kernels.cu

sptree.o: sptree.cpp
	$(NVCC)  $(INCLUDES) $(FLAGS) -c -o sptree.o sptree.cpp

t_sne_gpu.o: t_sne_gpu.cpp
	$(NVCC)  $(INCLUDES) $(FLAGS) -c -o t_sne_gpu.o t_sne_gpu.cpp

clean:
	rm -f *.a *.o 
