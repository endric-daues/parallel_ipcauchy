# Parallel IPCauchy Implementation - Final Project Submission for APMA E4302 Fall 2020 - Endric Daues

The code requires the installation of the GSL library. The code must be compiled and linked with reference to the GSL library as follows:

Circular Path:
```
mpicc -Wall -o main_mpi main_mpi.c helpers.c integrand.c -I/moto/home/ed2691/gsl/include -L/moto/home/ed2691/gsl/bin -std=c99 -Wall -lgsl -lgslcblas -lm
```

Shortest Path:
```
mpicc -Wall -o sp_main_mpi sp_main_mpi.c shortest_path.c helpers.c integrand.c -I/moto/home/ed2691/gsl/include -L/moto/home/ed2691/gsl/bin -std=c99 -Wall -lgsl -lgslcblas -lm
```

The results from running the code on various processes/nodes can be found in the specific directories.



