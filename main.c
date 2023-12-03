#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "functions.h"

int main(int argc, char *argv[]) {

    // Ensuring 2 arguments
    if (argc != 2) {
        fprintf(stderr, "Sorry, I expected 2 arguments. Please use the following format: %s <filename.mtx>\n", argv[0]);
        return 1;
    }

    // Convert mtx file to CSRMatrix
    CSRMatrix *csrMatrix = (CSRMatrix*)malloc(sizeof(CSRMatrix));
    if (csrMatrix == NULL) {
        fprintf(stderr, "Error allocating memory for CSRMatrix.\n");
        exit(EXIT_FAILURE);
    }
    ReadMMtoCSR(argv[1], csrMatrix);

    // Print the CSRMatrix
    // printCSRMatrix(csrMatrix);

    // Create a PNG image representing the sparsity pattern of the matrix
    createSparsePatternImage(csrMatrix->col_ind, csrMatrix->row_ptr, csrMatrix->num_rows, csrMatrix->num_cols, argv[1]);

    // Initialize b and x vectors
    double* b = (double*)malloc(csrMatrix->num_rows * sizeof(double)); // RHS vector initialization
    double* x = (double*)malloc(csrMatrix->num_rows * sizeof(double)); // Solution vector initialization
    if (b == NULL) {
        fprintf(stderr, "Error allocating memory for b vector.\n");
        exit(EXIT_FAILURE);
    } else if (x == NULL) {
        fprintf(stderr, "Error allocating memory for x vector.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < csrMatrix->num_rows; i++) {
        b[i] = 1.0; // Assume b = [1, 1, ...]
        x[i] = 0.0; // Initial BiCGSTAB guess is x = [0, 0, ...]
    }

    // Solve via BiCGSTAB and time it
    clock_t start = clock();
    bicgstab(csrMatrix, b, x, 1e-7, 10000);
    
    // Compute BiCGSTAB residual norm
    double* Ax = (double*)malloc(csrMatrix->num_rows * sizeof(double));
    spmv_csr(csrMatrix, x, Ax);
    double bicgstab_residual = 0.0;
    for (int i = 0; i < csrMatrix->num_rows; i++) {
        bicgstab_residual += (Ax[i] - b[i]) * (Ax[i] - b[i]); // Residual = Ax - b
        // printf("%lf ", x[i]); // Use this to print x
    }
    bicgstab_residual = sqrt(bicgstab_residual);

    // Solve via Conjugate Gradient
    conjugate_gradient(csrMatrix, b, x, 1e-7, 10000);
    Ax = (double*)realloc(Ax, csrMatrix->num_rows * sizeof(double));
    spmv_csr(csrMatrix, x, Ax);
    double conj_grad_residual = 0.0;
    for (int i = 0; i < csrMatrix->num_rows; i++) {
        conj_grad_residual += (Ax[i] - b[i]) * (Ax[i] - b[i]); // Residual = Ax - b
        // printf("%lf ", x[i]); // Use this to print x
    }
    conj_grad_residual = sqrt(conj_grad_residual);

    double optimal_residual = (bicgstab_residual < conj_grad_residual) ? bicgstab_residual : conj_grad_residual;

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    printf("\nMatrix Name: %s\n", argv[1]);
    printf("Matrix Dimensions: %d x %d\n", csrMatrix->num_rows, csrMatrix->num_cols);
    printf("Number of Non-Zero Values: %d\n", csrMatrix->num_non_zeros);
    printf("Program Runtime: %f seconds\n", time_spent);
    // printf("Biconjugate Gradient Stabilzed Residual Norm: %f\n", bicgstab_residual);
    // printf("Conjugate Gradient Residual Norm: %f\n", conj_grad_residual);
    printf("Residual Norm: %f\n", optimal_residual);

    free(b);
    free(x);
    free(Ax);
    free(csrMatrix->csr_data);
    free(csrMatrix->col_ind);
    free(csrMatrix->row_ptr);
    free(csrMatrix);

    return 0;
}