all: serial.out parallel.out

parallel.out: parallel_rFFT.c
      gcc -o parallel.out parallel_rFFT.c -lm -fopenmp

serial.out: serial_rFFT.c
      gcc -o serial.out serial_rFFT.c -lm
      
      
