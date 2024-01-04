# Biconjugate Gradient Stabilized Algorithm
This README.md file is a transcribed from LaTeX and summarized. To view the full PDF version, click [here](https://github.com/andrewderango/Biconjugate-Gradient-Stabilized-Algorithm/blob/main/Extra%20README%20Files/README.pdf).

## 1: Introduction

### 1.1: Assignment Overview
The goal of this assignment was to develop a program in C that is capable of solving first-order
linear systems of equations. These systems are of the form Ax=b, where A is a known matrix, b
is a known vector, and x is a vector whose elements must be solved for. The algorithm should be
capable to solving very large and sparse matrices, up to millions of rows. In developing a solution,
two parameters must be optimized/minimized: The runtime and the norm of the residual vector.
The residual norm, or rather the magnitude of the residual vector, was the measure used to
illustrate the accuracy of the solution given and can be defined by ||Ax − b||.\
\
It should be noted that very large matrices can impose significant memory issues when running
the program. As the matrices being dealt with are sparse, it is better to store the matrices in CSR
format, rather than storing every element of the 2D array in memory. This not only minimizes the
memory allocated and reduces the incidence of segmentation fault, but it also improves runtime.
Furthermore, it would take a very long time to read the given matrices if each of their elements
were given in a file. This would entail reading trillions of values, which is simply infeasible. To
overcome this issue, files were provided and read from MTX format. For further information
regarding MTX formatting, visit the Given Files folder of this repo.

### 1.2: Approach
The general process that the program developed for this assignment undergoes is described in this
section.\
\
First, it reads information from the MTX file and converts it directly to the CSR format.
Then, it checks if the matrix is triangular. If the matrix is triangular, then it asks the user if they
intended to solve the linear system where A is a simply the triangular matrix as provided, or if
A was intended to be the corresponding symmetrical matrix. This occurs because a triangular
matrix and its symmetrical counterpart, obtained by reflecting one triangular obliquely onto the
other triangle, would have the same MTX file. In MTX format, there is no way to discern the 
triangular matrix from the symmetrical. If the user decides to symmetrise the matrix, or reflect
the triangle obliquely to the other triangle, then the program completes the rest of the operations
using the symmetrical matrix.\
\
Next, the program uses the ```png.h``` library to generate a PNG file that shows the sparsity of
the matrix. This is an image that shows the sparsity pattern of the matrix, represented by black
and white pixels showing where non-zero values exist in the matrix. White represents non-zero
elements, while black pixels represent zeros. It creates a new directory if one has not already been
created to add the PNG file to. The PNG file is always added to a folder in the directory that
the user is currently in called Sparsity Pattern Images. The program notifies the user of this and
specifies the file name and directory.\
\
Following this, the program undergoes the iterative Biconjugate Gradient Stabilized (BiCGSTAB)
algorithm and returns a solution vector x. It then computes and stores the residual and its norm.
The program then executes the Conjugate Gradient algorithm and stores its residual norm. The
program then compares the residual norms and takes the better (lower) one. The program calculates the program runtime, which accounts for both algorithms combined.

### 1.3: Biconjugate Gradient Stabilized Algorithm Overview
BiCGSTAB is a robust algorithm used to solve first-order linear systems of equations. It is an
extension of the standard conjugate graduate method, in that it supports solutions for matrices
that are both nonsymmetrical and indefinite. That is, the eigenvalues of the A matrix are not
required to be positive in order to derive an accurate solution to the system. As an iterative
method that uses the residual to define two new search directions (biconjugate gradient vectors),
BiCGSTAB also has functionality that enforces stabilization which renders faster convergence.

## 2: Implementation in C

### 2.1: Structure Overview
The code can be broken down into four files:\
• ```Makefile```: Defines rules to be employed upon compilation.\
• ```functions.h```: Contains function prototypes that are formally defined in ```functions.c```\
• ```main.c```: Contains the main function that calls functions from ```functions.c```\
• ```functions.c```: Defines all functions other than int main such as ```spvm_csr```, ```bicgstab```,
```conjugate_gradient```, and more.

### 2.2: Parameter Adjustment
There exists four arbitrary parameters within the program that may significantly change the resultant solution for ```x```, and thus the residual vector and its norm. They are defined in the following
lines:
```
bicgstab ( csrMatrix , b , x , 1e -7 , 10000);
conjugate_gradient ( csrMatrix , b , x , 1e -7 , 10000);
```
Both the BiCGSTAB and Conjugate Gradient algorithms take in the parameters ```tolerance``` and
```max_iterations```. Above, they are defined as 0.0000001 and 10,000 respectively.\
• ```tolerance```: The iterations stop once the residual norm converges to below this specified
tolerance. This should depend on the accuracy that the user is looking for in their specific
circumstance.\
• ```max_iterations```: If the algorithms can’t converge, then the iterations will stop after max iterations
iterations. This helps deal with cases in which BiCGSTAB will take a long time to converge
to the true solution.

### 2.3: CSR Formatting
The CSR format of LFAT5.mtx can be seen below:
```
Number of Non-Zeros: 46
Row Pointer: 0 3 5 7 11 15 18 21 26 31 33 35 39 43 46
Column Index: 0 3 4 1 5 2 6 0 3 7 8 0 4 7 8 1 5 9 2 6 10 3 4 7 11 12 3 4 8 11 12 5 9 6 10 7 8 11 13 7 8 12 13 11 12 13
CSR Data: 1.570880 -94.252800 0.785440 12566400.000000 -6283200.000000 0.608806 -0.304403 -94.252800 15080.448000 -7540.224000 94.252800 0.785440 3.141760 -94.252800 0.785440 -6283200.000000 12566400.000000 -6283200.000000 -0.304403 0.608806 -0.304403 -7540.224000 -94.252800 15080.448000 -7540.224000 94.252800 94.252800 0.785440 3.141760 -94.252800 0.785440 -6283200.000000 12566400.000000 -0.304403 0.608806 -7540.224000 -94.252800 15080.448000 94.252800 94.252800 0.785440 3.141760 0.785440 94.252800 0.785440 1.570880
```

### 2.4: Dependencies
This program requires ```png.h``` in order to execute.
\
To install the dependency on Debian-based Linux, run the following:
```
sudo apt-get update
sudo apt-get install libpng-dev
```
On macOS, use Homebrew:
```
brew update
brew install libpng
```

### 2.5: Running the Program
To run the program, the following commands must be run in the program’s directory:
```
make
./bicgstab <filename.mtx>
```
For example, if you want to solve the system Ax=b for b = [1, 1, ..., 1] and A = [LFAT5.mtx](https://github.com/andrewderango/Biconjugate-Gradient-Stabilized-Algorithm/blob/main/Matrices/LFAT5.mtx), then enter the following:
```
make
./bicgstab LFAT5.mtx
```

### 2.6: Sample Output
A sample terminal output for LFAT5.mtx can be seen below:
```
andrewderango@Andrews-MacBook-Air-2 Biconjugate-Gradient-Stabilized-Algorithm % ./bicgstab LFAT5.mtx
I have detected that this matrix is upper triangular. Would you like to symmetrise it? (Y/N): y
Sparsity pattern image saved to your current directory: Sparsity Pattern Images/LFAT5_sparsity_pattern.png

Matrix Name: LFAT5.mtx
Matrix Dimensions: 14 x 14
Number of Non-Zero Values: 46
Program Runtime: 0.000060 seconds
Residual Norm: 0.000000
```

## 3: Results

| Matrix              | Dimensions       | Non-Zeros  | CPU Time (s) | Residual Norm |
|:-------------------:|:----------------:|:----------:|:------------:|:-------------:|
| b1_ss.mtx          | 7x7              | 15         | 0.000010     | 0.000000      |
| LFAT5.mtx           | 14x14            | 46         | 0.000038     | 0.000000      |
| LF10.mtx            | 18x18            | 82         | 0.000141     | 0.000000      |
| ex3.mtx             | 1821x1821        | 52685      | 1.011970     | 0.000112      |
| jnlbrng1.mtx        | 40000x40000      | 199200     | 0.041605     | 0.000000      |
