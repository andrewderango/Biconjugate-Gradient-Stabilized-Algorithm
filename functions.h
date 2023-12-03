#ifndef FUNCTIONS_H
#define FUNCTIONS_H

typedef struct {
    double *csr_data;   // Array of non-zero values
    int *col_ind;       // Array of column indices
    int *row_ptr;       // Array of row pointers
    int num_non_zeros;  // Number of non-zero elements
    int num_rows;       // Number of rows in matrix
    int num_cols;       // Number of columns in matrix
} CSRMatrix;

typedef struct {
    int row;            // Index of row in matrix
    int col;            // Index of column in matrix
    double value;       // Value of specific element in matrix
} MTXRow;

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix);
int compareMTXData(const void* a, const void* b);
void createSparsePatternImage(const int *columns, const int *row_ptr, int num_rows, int num_cols, const char *mtx_file_name);
void spmv_csr(const CSRMatrix *A, const double *x, double *y);
void bicgstab(CSRMatrix *A, double *b, double *x, double tolerance, int max_iterations);
void conjugate_gradient(CSRMatrix *A, double *b, double *x, double tolerance, int max_iterations);

#endif
