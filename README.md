# MPI_wrapper

This is a basic 2D MPI wrapper for integration into any pre-existing Fortran code. It contains the beginnings of a much more complicated MPI scheme, but for the moment it is assumed that no communication between CPUs is required, as all tasks are treated as independent. Adding extra dimensions is very simple to create a 3D or 4D MPI scheme, if necessary.

Work on this basic scheme will allow one to develop their skills in parallel computing as they learn how to send and receive data between CPUs, discover how to map a multi-CPU task and how each CPU will interpret the source code, and collate output data on the master CPU, or integrate MPI/IO into it.

I hope you find it a useful start to your path into parallel computing with Fortran.