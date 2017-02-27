all: dgemm_cpu dgemm_gpu

dgemm_cpu: dgemm_test.F90 timer.o
	pgf90 dgemm_test.F90 timer.o -o dgemm_cpu -lblas

dgemm_gpu: dgemm_test.F90 timer.o
	pgf90 dgemm_test.F90 timer.o -o dgemm_gpu -lblas -Mcuda=cc35,cuda8.0 -D__CUDA__

timer.o: timer.c
	gcc -c timer.c

clean:
	rm -f *.o dgemm_gpu dgemm_cpu 
