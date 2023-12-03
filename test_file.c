#include <stdio.h>
#include <stdlib.h>
// #include <stdbool.h>
#include <math.h>
#include <time.h>

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

void spmv_csr(const CSRMatrix *A, const double *x, double *y) {
    for (int i = 0; i < A->num_rows; i++) {
        y[i] = 0.0;
        for (int j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++) {
            y[i] += A->csr_data[j] * x[A->col_ind[j]];
        }
    }
}

void printCSRMatrix(CSRMatrix* matrix) {
    printf("CSR Matrix:\n");
    printf("Number of Rows: %d\n", matrix->num_rows);
    printf("Number of Columns: %d\n", matrix->num_cols);
    printf("Number of Non-Zero Elements: %d\n", matrix->num_non_zeros);
    
    // printf("CSR Data: ");
    // for (int i = 0; i < matrix->num_non_zeros; i++) {
    //     printf("%lf ", matrix->csr_data[i]);
    // }
    // printf("\n");
    
    // printf("Column Indices: ");
    // for (int i = 0; i < matrix->num_non_zeros; i++) {
    //     printf("%d ", matrix->col_ind[i]);
    // }
    // printf("\n");
    
    // printf("Row Pointers: ");
    // for (int i = 0; i < matrix->num_rows + 1; i++) {
    //     printf("%d ", matrix->row_ptr[i]);
    // }
    // printf("\n");

    printf("Matrix:\n");
    for (int i = 0; i < matrix->num_rows; i++) {
        for (int j = 0; j < matrix->num_cols; j++) {
            int index = -1;
            for (int k = matrix->row_ptr[i]; k < matrix->row_ptr[i+1]; k++) {
                if (matrix->col_ind[k] == j) {
                    index = k;
                    break;
                }
            }
            if (index != -1) {
                printf("%.2lf\t", matrix->csr_data[index]);
            } else {
                printf("0.0\t");
            }
        }
        printf("\n");
    }
}

// Comparison function for MTXRow
int compareMTXData(const void* a, const void* b) {
    MTXRow* rowA = (MTXRow*)a;
    MTXRow* rowB = (MTXRow*)b;

    // Compare rows
    if (rowA->row != rowB->row) {
        return rowA->row - rowB->row;
    }

    // Rows are equal, so compare columns
    return rowA->col - rowB->col;
}

void bicgstab(CSRMatrix* A, double* b, double* x, double tolerance, int max_iterations) {
    double* Ax = (double*)malloc(A->num_rows * sizeof(double)); // LHS vector
    double* r = (double*)malloc(A->num_rows * sizeof(double)); // Residual vector
    double* r_hat = (double*)malloc(A->num_rows * sizeof(double)); // Biorthogonalized residual vector
    double rho = 1.0; // Dot product of residual and biorthogonalized residual
    double alpha = 1.0; // Step size for search direction. Alpha = rho/(rhat*v)
    double omega = 1.0; // Step size for stabilization. Omega = (rhat*r)/(rhat*v)
    double* v = (double*)malloc(A->num_rows * sizeof(double)); // Temporary storage of A * p
    double* p = (double*)malloc(A->num_rows * sizeof(double)); // Search direction vector

    spmv_csr(A, x, Ax);

    printf("\nAx: ");
    for (int i = 0; i < A->num_rows; i++) {
        printf("%lf ", Ax[i]);
    }
    // for (int i = 0; i < A->num_rows; i++) {
    //     r[i] = b[i] - Ax[i];
    //     r_hat[i] = r[i];
    //     p[i] = r[i];
    // }

    printf("\nb: ");
    for (int i = 0; i < A->num_rows; i++) {
        printf("%lf ", b[i]);
    }
    printf("\nr: ");
    for (int i = 0; i < A->num_rows; i++) {
        r[i] = b[i] - Ax[i];
        r_hat[i] = r[i];
        v[i] = 0.0;
        p[i] = 0.0;
        printf("%lf ", r[i]);
    }
    printf("\nn: %d\n", A->num_rows);

    printf("r_hat: ");
    for (int i = 0; i < A->num_rows; i++) {
        printf("%lf ", r_hat[i]);
    }
    printf("\n\n");

    for (int iteration = 0; iteration < max_iterations; iteration++) {
        printf("-- ITERATION %d --\n", iteration+1);

        double rho_new = 0.0;
        for (int i = 0; i < A->num_rows; i++) {
            rho_new += r_hat[i] * r[i];
        }
        printf("rho_new: %f\n", rho_new);

        double beta = (rho_new / (rho + 1e-16)) * (alpha / (omega + 1e-16));
        rho = rho_new;
        printf("beta: %f\n", beta);

        for (int i = 0; i < A->num_rows; i++) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        printf("p: ");
        for (int i = 0; i < A->num_rows; i++) {
            printf("%f ", p[i]);
        }
        printf("\n");

        spmv_csr(A, p, v);
        printf("v: ");
        for (int i = 0; i < A->num_rows; i++) {
            printf("%f ", v[i]);
        }
        printf("\n");

        double dot_r_hat_v = 0.0;
        for (int i = 0; i < A->num_rows; i++) {
            dot_r_hat_v += r_hat[i] * v[i];
        }
        alpha = rho / dot_r_hat_v;
        printf("alpha: %f\n", alpha);

        double* s = (double*)malloc(A->num_rows * sizeof(double)); // Temporary storage r - alpha * v
        double* t = (double*)malloc(A->num_rows * sizeof(double)); // Temporary storage A * s
        for (int i = 0; i < A->num_rows; i++) {
            s[i] = r[i] - alpha * v[i];
        }
        printf("s: ");
        for (int i = 0; i < A->num_rows; i++) {
            printf("%f ", s[i]);
        }
        printf("\n");

        spmv_csr(A, s, t);
        printf("t: ");
        for (int i = 0; i < A->num_rows; i++) {
            printf("%f ", t[i]);
        }
        printf("\n");

        double dot_t_s = 0.0;
        double dot_t_t = 0.0;
        for (int i = 0; i < A->num_rows; i++) {
            dot_t_s += t[i] * s[i];
            dot_t_t += t[i] * t[i];
        }

        omega = dot_t_s / (dot_t_t + 1e-16);
        printf("omega: %f\n", omega);

        for (int i = 0; i < A->num_rows; i++) {
            x[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }

        printf("x: ");
        for (int i = 0; i < A->num_rows; i++) {
            printf("%f ", x[i]);
        }
        printf("\n");

        printf("r: ");
        for (int i = 0; i < A->num_rows; i++) {
            printf("%f ", r[i]);
        }
        printf("\n");

        // Calculate norm of r
        double residual = 0.0;
        for (int i = 0; i < A->num_rows; i++) {
            residual += r[i] * r[i];
        }
        residual = sqrt(residual);
        printf("Residual: %f\n", residual);

        if (residual < tolerance) {
            break;
        }

        free(s);
        free(t);
    }

    free(Ax);
    free(r);
    free(r_hat);
    free(v);
    free(p);
}

CSRMatrix* convert_to_csr_from_mtx(const char* filename) {
    // Try to open file
    FILE *file = NULL;
    file = fopen(filename, "r"); // Returns pointer to file if successful, NULL otherwise

    // If the file doesn't exist, return NULL
    if (file == NULL) {
        return NULL;
    }

    int rows, columns, nonzero_values, *row_ind;
    double mtx_row, mtx_column, mtx_value;
    int lower_triangular = 1, upper_triangular = 1;

    // Skip all the headers starting the line with %
    char line[1024]; // Will fail if a line in mtx > 1024 chars ###
    line[0] = '%';
    while (line[0] == '%') {
        fgets(line, sizeof(line), file);
    };
    sscanf(line, "%d %d %d", &rows, &columns, &nonzero_values); // First line of mtx has rows, columns, nonzero values

    // Allocate memory for the CSRMatrix
    CSRMatrix* matrix = (CSRMatrix*)malloc(sizeof(CSRMatrix));
    matrix->num_rows = rows;
    matrix->num_cols = columns;
    matrix->num_non_zeros = nonzero_values;

    matrix->csr_data = (double*)malloc(nonzero_values * sizeof(double));
    matrix->col_ind = (int*)malloc(nonzero_values * sizeof(int));
    matrix->row_ptr = (int*)malloc((rows + 1) * sizeof(int));
    row_ind = (int*)malloc(nonzero_values * sizeof(int));

    // Initialize row_ptr to 0
    for (int i = 0; i < rows+1; i++) {
        matrix->row_ptr[i] = 0;
    }

    // Read the values into the CSRMatrix
    int current_row = 0;
    int max_col = 0;
    for (int i = 0; i < matrix->num_non_zeros; i++) {
        fscanf(file, "%lf %lf %lf", &mtx_row, &mtx_column, &mtx_value);
        // printf("%lf %lf %lf\n", mtx_row, mtx_column, mtx_value);

        if (mtx_column < max_col) {
            printf("This is improper MTX format. MTX format requires the column-index column to be ascending.\n");
            return NULL;
        } else {
            max_col = mtx_column;
        }

        if (mtx_value != 0.0) {
            matrix->csr_data[i] = mtx_value;
            matrix->col_ind[i] = mtx_column - 1;
            row_ind[i] = mtx_row - 1;

            if (mtx_row > mtx_column) {
                lower_triangular = 0;
            } else if (mtx_row < mtx_column) {
                upper_triangular = 0;
            }
        }
    }

    // Initialize row_ptr to 0
    for (int i = 0; i < rows+1; i++) {
        matrix->row_ptr[i] = 0;
    }

    // Count the number of non-zero elements in each row
    for (int i = 0; i < matrix->num_non_zeros; i++) {
        matrix->row_ptr[row_ind[i] + 1]++;
    }

    // Compute the cumulative sum of row_ptr
    for (int i = 0; i < rows; i++) {
        matrix->row_ptr[i + 1] += matrix->row_ptr[i];
    }

    MTXRow* data = (MTXRow*)malloc(matrix->num_non_zeros * sizeof(MTXRow));
    for (int i = 0; i < matrix->num_non_zeros; i++) {
        data[i].row = row_ind[i];
        data[i].col = matrix->col_ind[i];
        data[i].value = matrix->csr_data[i];
    }

    qsort(data, matrix->num_non_zeros, sizeof(MTXRow), compareMTXData);

    // Extract the sorted data
    for (int i = 0; i < matrix->num_non_zeros; i++) {
        row_ind[i] = data[i].row;
        matrix->col_ind[i] = data[i].col;
        matrix->csr_data[i] = data[i].value;
    }

    // Free the array
    free(data);

    fclose(file);

    // printf("\n");
    // printf("Num Rows: %d\n", matrix->num_rows);
    // printf("Num Cols: %d\n", matrix->num_cols);
    // printf("Num Non-Zeros: %d\n", matrix->num_non_zeros);

    // printf("CSR Data: ");
    // for (int i = 0; i < matrix->num_non_zeros; i++) {
    //     printf("%lf ", matrix->csr_data[i]);
    // }
    // printf("\n");
    
    // printf("Column Indices: ");
    // for (int i = 0; i < matrix->num_non_zeros; i++) {
    //     printf("%d ", matrix->col_ind[i]);
    // }
    // printf("\n");

    // printf("Row Indices: ");
    // for (int i = 0; i < matrix->num_non_zeros; i++) {
    //     printf("%d ", row_ind[i]);
    // }
    // printf("\n");
    
    // printf("Row Pointers: ");
    // for (int i = 0; i < matrix->num_rows + 1; i++) {
    //     printf("%d ", matrix->row_ptr[i]);
    // }
    // printf("\n");

    // Ask user if they want to symmetrise the matrix if it is triangular
    if (upper_triangular == 1 && lower_triangular == 0) {
        printf("I have detected that this matrix is upper triangular. Would you like to symmetrise it? (Y/N): ");
    } else if (upper_triangular == 0 && lower_triangular == 1) {
        printf("I have detected that this matrix is lower triangular. Would you like to symmetrise it? (Y/N): ");
    }

    if (upper_triangular == 1 || lower_triangular == 1) {
        char symmetrise;
        scanf(" %c", &symmetrise);
        while (symmetrise != 'y' && symmetrise != 'Y' && symmetrise != '1' && symmetrise != 'n' && symmetrise != 'N' && symmetrise != '0') {
            printf("Invalid response. Please enter Y or N: ");
            scanf(" %c", &symmetrise);
        }
        if (symmetrise == 'y' || symmetrise == 'Y' || symmetrise == '1') {
            int og_num_non_zeros = matrix->num_non_zeros;
            for (int i = 0; i < og_num_non_zeros; i++) {
                if (matrix->col_ind[i] != row_ind[i]) {
                    matrix->num_non_zeros++;
                    matrix->col_ind = realloc(matrix->col_ind, matrix->num_non_zeros * sizeof(int));
                    row_ind = realloc(row_ind, matrix->num_non_zeros * sizeof(int));
                    matrix->csr_data = realloc(matrix->csr_data, matrix->num_non_zeros * sizeof(double));
                    matrix->col_ind[matrix->num_non_zeros - 1] = row_ind[i];
                    row_ind[matrix->num_non_zeros - 1] = matrix->col_ind[i];
                    matrix->csr_data[matrix->num_non_zeros - 1] = matrix->csr_data[i];
                }
            }

            // Re-sort CSR data
            MTXRow* data = (MTXRow*)malloc(matrix->num_non_zeros * sizeof(MTXRow));
            for (int i = 0; i < matrix->num_non_zeros; i++) {
                data[i].row = row_ind[i];
                data[i].col = matrix->col_ind[i];
                data[i].value = matrix->csr_data[i];
            }
            qsort(data, matrix->num_non_zeros, sizeof(MTXRow), compareMTXData);
            // Extract the sorted data
            for (int i = 0; i < matrix->num_non_zeros; i++) {
                row_ind[i] = data[i].row;
                matrix->col_ind[i] = data[i].col;
                matrix->csr_data[i] = data[i].value;
            }
            free(data);
            // Re-compute row_ptr
            for (int i = 0; i < rows+1; i++) {
                matrix->row_ptr[i] = 0;
            }
            for (int i = 0; i < matrix->num_non_zeros; i++) {
                matrix->row_ptr[row_ind[i] + 1]++;
            }
            for (int i = 0; i < rows; i++) {
                matrix->row_ptr[i + 1] += matrix->row_ptr[i];
            }
        }
    }

    printf("CSR Data: ");
    for (int i = 0; i < matrix->num_non_zeros; i++) {
        printf("%lf ", matrix->csr_data[i]);
    }
    printf("\n");
    
    printf("Column Indices: ");
    for (int i = 0; i < matrix->num_non_zeros; i++) {
        printf("%d ", matrix->col_ind[i]);
    }
    printf("\n");

    printf("Row Indices: ");
    for (int i = 0; i < matrix->num_non_zeros; i++) {
        printf("%d ", row_ind[i]);
    }
    printf("\n");

    printf("Row Pointers: ");
    for (int i = 0; i < matrix->num_rows + 1; i++) {
        printf("%d ", matrix->row_ptr[i]);
    }
    printf("\n");

    free(row_ind);

    return matrix;
}

int main(int argc, char *argv[]) {

    // Ensuring 2 arguments
    if (argc != 2) {
        fprintf(stderr, "Sorry, I expected 2 arguments. Please use the following format: %s <filename.mtx>\n", argv[0]);
        return 1;
    }

    // Convert mtx file to CSRMatrix
    CSRMatrix* csrMatrix = convert_to_csr_from_mtx(argv[1]);

    // If the conversion failed, return an error
    if (csrMatrix == NULL) {
        fprintf(stderr, "Sorry, the file %s could not be opened or converted to the CSR format.\n", argv[1]);
        return 1;
    }

    // Print the CSRMatrix
    // printCSRMatrix(csrMatrix);

    // Assume b is [1, 1, ..., 1]
    double* b = (double*)malloc(csrMatrix->num_rows * sizeof(double));
    double* x = (double*)malloc(csrMatrix->num_rows * sizeof(double)); // Solution vector initialization
    for (int i = 0; i < csrMatrix->num_rows; i++) {
        b[i] = 1.0;
        x[i] = 1.0;
    }

    // Solve via BiCGSTAB
    bicgstab(csrMatrix, b, x, 1e-6, 1000);

    // Free the memory
    free(b);
    free(x);
    free(csrMatrix->csr_data);
    free(csrMatrix->col_ind);
    free(csrMatrix->row_ptr);
    free(csrMatrix);

    return 0;
}

// TODO:
// check ###
// late: cleanup