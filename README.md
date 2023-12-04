# 1 Introduction
## 1.1 Assignment Overview
The goal of this assignment was to develop a program in C that is capable of solving first-order
linear systems of equations. These systems are of the form Ax=b, where A is a known matrix, b
is a known vector, and x is a vector whose elements must be solved for. The algorithm should be
capable to solving very large and sparse matrices, up to millions of rows. In developing a solution,
two parameters must be optimized/minimized: The runtime and the norm of the residual vector.
The residual norm, or rather the magnitude of the residual vector, was the measure used to
illustrate the accuracy of the solution given and can be defined by ||Ax − b||.
It should be noted that very large matrices can impose significant memory issues when running
the program. As the matrices being dealt with are sparse, it is better to store the matrices in CSR
format, rather than storing every element of the 2D array in memory. This not only minimizes the
memory allocated and reduces the incidence of segmentation fault, but it also improves runtime.
Furthermore, it would take a very long time to read the given matrices if each of their elements
were given in a file. This would entail reading trillions of values, which is simply infeasible. To
overcome this issue, files were provided and read from MTX format. For further information
regarding MTX formatting, visit the Given Files folder of this repo.
