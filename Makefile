all:
	nvcc -arch=sm_52 -O3 -dc *.cu
	nvcc -arch=sm_52 -O3 *.o -o pfunc-cuda

clean:
	rm -f pfunc-cuda
	rm -f *.o

debug:
	nvcc -arch=sm_52 -G -dc *.cu
	nvcc -arch=sm_52 -G *.o -o pfunc-cuda

