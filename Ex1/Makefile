all: vadd_cpu vadd_gpu

vadd_cpu: main.F90 vector_add.F90
	pgf90  main.F90 vector_add.F90 -o vadd_cpu

vadd_gpu: main.F90 vector_add.F90
	
	pgf90 main.F90 vector_add.F90 -o vadd_gpu -lblas -Mcuda=cc35,cuda8.0 -D__CUDA__

clean:
	rm *.o *.mod vadd_cpu vadd_gpu
