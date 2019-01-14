all:
	nvcc -arch=sm_52 -O3 -dc -Xcompiler -fPIC *.cu
	nvcc -arch=sm_52 --shared -o cupyck.so *.o

clean:
	rm -f cupyck.so
	rm -f *.o

debug:
	nvcc -arch=sm_52 -g -G -dc -Xcompiler -fPIC *.cu
	nvcc -arch=sm_52 -g -G --shared -o cupyck.so *.o
