#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <png.h>
#include <string.h>

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

void printCSRMatrix(CSRMatrix *matrix) {
    printf("CSR Matrix:\n");
    printf("Number of Rows: %d\n", matrix->num_rows);
    printf("Number of Columns: %d\n", matrix->num_cols);
    printf("Number of Non-Zero Elements: %d\n", matrix->num_non_zeros);

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
int compareMTXData(const void *a, const void *b) {
    MTXRow* rowA = (MTXRow*)a;
    MTXRow* rowB = (MTXRow*)b;

    // Compare rows
    if (rowA->row != rowB->row) {
        return rowA->row - rowB->row;
    } else {
        return rowA->col - rowB->col; // compare columns if rows are equal
    }
}

void bicgstab(CSRMatrix *A, double *b, double *x, double tolerance, int max_iterations) {
    double* Ax = (double*)malloc(A->num_rows * sizeof(double)); // LHS vector
    double* r = (double*)malloc(A->num_rows * sizeof(double)); // Residual vector
    double* r_hat = (double*)malloc(A->num_rows * sizeof(double)); // Biorthogonalized residual vector
    double rho = 1.0; // Dot product of residual and biorthogonalized residual
    double alpha = 1.0; // Step size for search direction. Alpha = rho/(rhat*v)
    double omega = 1.0; // Step size for stabilization. Omega = (rhat*r)/(rhat*v)
    double* v = (double*)malloc(A->num_rows * sizeof(double)); // Temporary storage of A * p
    double* p = (double*)malloc(A->num_rows * sizeof(double)); // Search direction vector
    double best_residual; // Best residual so far
    double* best_x = (double*)malloc(A->num_rows * sizeof(double)); // Best solution so far

    if (Ax == NULL || r == NULL || r_hat == NULL || v == NULL || p == NULL || best_x == NULL) {
        fprintf(stderr, "Error allocating memory for vectors in BiCGSTAB function.\n");
        exit(EXIT_FAILURE);
    }

    spmv_csr(A, x, Ax);

    for (int i = 0; i < A->num_rows; i++) {
        r[i] = b[i] - Ax[i];
        r_hat[i] = r[i];
        v[i] = 0.0;
        p[i] = 0.0;
    }

    for (int iteration = 0; iteration < max_iterations; iteration++) {
        double rho_new = 0.0;
        for (int i = 0; i < A->num_rows; i++) {
            rho_new += r_hat[i] * r[i];
        }
        double beta = (rho_new / (rho + 1e-16)) * (alpha / (omega + 1e-16));
        rho = rho_new;
        for (int i = 0; i < A->num_rows; i++) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        spmv_csr(A, p, v);
        double dot_r_hat_v = 0.0;
        for (int i = 0; i < A->num_rows; i++) {
            dot_r_hat_v += r_hat[i] * v[i];
        }
        alpha = rho / dot_r_hat_v;
        double* s = (double*)malloc(A->num_rows * sizeof(double)); // Temporary storage r - alpha * v
        double* t = (double*)malloc(A->num_rows * sizeof(double)); // Temporary storage A * s
        if (s == NULL || t == NULL) {
            fprintf(stderr, "Error allocating memory for s or t vector.\n");
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < A->num_rows; i++) {
            s[i] = r[i] - alpha * v[i];
        }
        spmv_csr(A, s, t);
        double dot_t_s = 0.0;
        double dot_t_t = 0.0;
        for (int i = 0; i < A->num_rows; i++) {
            dot_t_s += t[i] * s[i];
            dot_t_t += t[i] * t[i];
        }
        omega = dot_t_s / (dot_t_t + 1e-16);
        for (int i = 0; i < A->num_rows; i++) {
            x[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }

        double residual = 0.0;
        for (int i = 0; i < A->num_rows; i++) {
            residual += r[i] * r[i];
        }
        residual = sqrt(residual);
        if (iteration == 0 || residual < best_residual) {
            best_residual = residual;
            for (int i = 0; i < A->num_rows; i++) {
                best_x[i] = x[i];
            }
        }

        // printf("Iteration: %d\tResidual: %lf\n", iteration, residual);

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

    for (int i = 0; i < A->num_rows; i++) {
        x[i] = best_x[i];
    }
    // printf("Residual: %f\n", best_residual);
}

void conjugate_gradient(CSRMatrix* A, double* b, double* x, double tolerance, int max_iterations) {
    double* r = (double*)malloc(A->num_rows * sizeof(double)); // Residual vector
    double* p = (double*)malloc(A->num_rows * sizeof(double)); // Search direction vector
    double* Ap = (double*)malloc(A->num_rows * sizeof(double)); // Temporary storage of A * p
    double alpha, beta;
    double best_residual;
    double* best_x = (double*)malloc(A->num_rows * sizeof(double));

    if (r == NULL || p == NULL || Ap == NULL || best_x == NULL) {
        fprintf(stderr, "Error allocating memory for vectors in conjugate gradient function.\n");
        exit(EXIT_FAILURE);
    }

    spmv_csr(A, x, Ap);

    for (int i = 0; i < A->num_rows; i++) {
        r[i] = b[i] - Ap[i];
        p[i] = r[i];
    }
    double r_dot_r = 0;
    for (int i = 0; i < A->num_rows; i++) {
        r_dot_r += r[i]*r[i];
    }

    for (int iteration = 0; iteration < max_iterations; iteration++) {
        spmv_csr(A, p, Ap);
        
        double p_dot_Ap = 0.0;
        for (int i = 0; i < A->num_rows; i++) {
            p_dot_Ap += p[i] * Ap[i];
        }
        double alpha = r_dot_r / (p_dot_Ap + 1e-16);
        for (int i = 0; i < A->num_rows; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        double r_dot_r_new = 0;
        for (int i = 0; i < A->num_rows; i++) {
            r_dot_r_new += r[i]*r[i];
        }
        double beta = r_dot_r_new / (r_dot_r + 1e-16);
        for (int i = 0; i < A->num_rows; i++) {
            p[i] = r[i] + beta * p[i];
        }
        r_dot_r = r_dot_r_new;
        double residual = 0.0;
        for (int i = 0; i < A->num_rows; i++) {
            residual += r[i] * r[i];
        }
        residual = sqrt(residual);

        if (iteration == 0 || residual < best_residual) {
            best_residual = residual;
            for (int i = 0; i < A->num_rows; i++) {
                best_x[i] = x[i];
            }
        }

        // printf("Iteration: %d\tResidual: %lf\n", iteration, residual);

        if (residual < tolerance) {
            residual = best_residual;
            break;
        }
    }

    free(r);
    free(p);
    free(Ap);

    for (int i = 0; i < A->num_rows; i++) {
        x[i] = best_x[i];
    }
}

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix) {
    // Try to open file
    FILE *file = NULL;
    file = fopen(filename, "r"); // Returns pointer to file if successful or NULL if unsuccessful

    // If the file doesn't exist or can't open it, print error message and quit the program
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
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

    matrix->num_rows = rows;
    matrix->num_cols = columns;
    matrix->num_non_zeros = nonzero_values;

    matrix->csr_data = (double*)malloc(nonzero_values * sizeof(double));
    matrix->col_ind = (int*)malloc(nonzero_values * sizeof(int));
    matrix->row_ptr = (int*)malloc((rows + 1) * sizeof(int));
    row_ind = (int*)malloc(nonzero_values * sizeof(int));

    if (matrix->csr_data == NULL || matrix->col_ind == NULL || matrix->row_ptr == NULL || row_ind == NULL) {
        fprintf(stderr, "Error allocating memory for CSRMatrix.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize row_ptr to 0
    for (int i = 0; i < rows+1; i++) {
        matrix->row_ptr[i] = 0;
    }

    // Read the values into the CSRMatrix
    int current_row = 0;
    int max_col = 0;
    for (int i = 0; i < matrix->num_non_zeros; i++) {
        fscanf(file, "%lf %lf %lf", &mtx_row, &mtx_column, &mtx_value);

        if (mtx_column < max_col) {
            printf("This is improper MTX format. MTX format requires the column-index column to be ascending.\n");
            exit(EXIT_FAILURE);
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
    if (data == NULL) {
        fprintf(stderr, "Error allocating memory for data.\n");
        exit(EXIT_FAILURE);
    }
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
    fclose(file);

    // Ask user if they want to symmetrise the matrix if it is triangular
    if (upper_triangular == 1 && lower_triangular == 0) {
        printf("I have detected that this matrix is upper triangular. Would you like to symmetrise it? (Y/N): ");
    } else if (upper_triangular == 0 && lower_triangular == 1) {
        printf("I have detected that this matrix is lower triangular. Would you like to symmetrise it? (Y/N): ");
    }

    if (upper_triangular == 1 || lower_triangular == 1) {
        char symmetrise;
        scanf(" %c", &symmetrise);
        while (symmetrise != 'y' && symmetrise != 'Y' && symmetrise != '1' && symmetrise != 'n' && symmetrise != 'N' && symmetrise != '0') { // fix ###
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
            if (data == NULL) {
                fprintf(stderr, "Error allocating memory for data.\n");
                exit(EXIT_FAILURE);
            }
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

    free(row_ind);
}

// Function to create a PNG image representing the sparsity pattern of a CSR matrix ###
void createSparsePatternImage(const double *values, const int *columns, const int *row_ptr, int num_rows, int num_cols, const char *mtx_file_name) {
    if (num_rows > 50000 || num_cols > 50000) {
        fprintf(stderr, "Sorry, the matrix is too large to create a sparsity pattern image.\n");
        return;
    }

    // Remove extension of input matrix file name
    char* dot = strchr(mtx_file_name, '.');
    size_t base_name_length = dot != NULL ? dot - mtx_file_name : strlen(mtx_file_name);
    char base_name[256];  // Adjust the size as needed
    strncpy(base_name, mtx_file_name, base_name_length);
    base_name[base_name_length] = '\0';  // Null-terminate the string
    
    // Image size and scale factor
    int scale_factor;
    if (num_rows * num_cols > 5000) {
        scale_factor = 1;
    } else if (num_rows * num_cols > 500) {
        scale_factor = 10;
    } else {
        scale_factor = 100;
    }
    uint64_t image_width = num_cols * scale_factor;
    uint64_t image_height = num_rows * scale_factor;

    // Create a buffer to hold pixel data
    uint8_t* pixel_data = (uint8_t*)calloc(image_width * image_height, sizeof(uint8_t));
    if (pixel_data == NULL) {
        fprintf(stderr, "Error allocating memory for pixel data.\n");
        exit(EXIT_FAILURE);
    }
    // Populate pixel_data based on the sparse matrix
    for (int i = 0; i < num_rows; ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            int col = columns[j];
            for (uint32_t di = 0; di < scale_factor; ++di) {
                for (uint32_t dj = 0; dj < scale_factor; ++dj) {
                    pixel_data[(i * scale_factor + di) * image_width + col * scale_factor + dj] = 255;  // Set pixel to white
                }
            }
        }
    }

    // Create an array of pointers to the rows in the image
    uint8_t** rows = (uint8_t**)malloc(image_height * sizeof(uint8_t*));
    if (rows == NULL) {
        fprintf(stderr, "Memory allocation error.\n");
        free(pixel_data);
        exit(EXIT_FAILURE);
    }
    for (uint32_t i = 0; i < image_height; ++i) {
        rows[i] = &pixel_data[i * image_width];
    }

    // Write the image data to a file
    char file_name[256];  // Adjust the size as needed
    sprintf(file_name, "%s_sparsity_pattern.png", base_name);
    FILE* fp = fopen(file_name, "wb");
    if (!fp) abort();

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) abort();

    png_infop info = png_create_info_struct(png);
    if (!info) abort();

    if (setjmp(png_jmpbuf(png))) abort();

    png_init_io(png, fp);

    // Write header
    png_set_IHDR(
        png,
        info,
        image_width, image_height,
        8,
        PNG_COLOR_TYPE_GRAY,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png, info);  // Use 'info' instead of 'rows'

    // Write image data
    png_write_image(png, rows);  // Use 'rows' here

    // Finish writing
    png_write_end(png, NULL);
    if (png && info)
        png_destroy_write_struct(&png, &info);
    if (fp)
        fclose(fp);

    // Free allocated memory
    free(rows);
    free(pixel_data);

    printf("Sparsity pattern image saved to your current directory: %s\n", file_name);
}

int main(int argc, char *argv[]) {

    // Ensuring 2 arguments
    if (argc != 2) {
        fprintf(stderr, "Sorry, I expected 2 arguments. Please use the following format: %s <filename.mtx>\n", argv[0]);
        return 1;
    }

    // Convert mtx file to CSRMatrix
    CSRMatrix* csrMatrix = (CSRMatrix*)malloc(sizeof(CSRMatrix));
    if (csrMatrix == NULL) {
        fprintf(stderr, "Error allocating memory for CSRMatrix.\n");
        exit(EXIT_FAILURE);
    }
    ReadMMtoCSR(argv[1], csrMatrix);

    // If the conversion failed, return an error
    if (csrMatrix == NULL) {
        fprintf(stderr, "Sorry, the file %s could not be opened or converted to the CSR format.\n", argv[1]);
        return 1;
    }

    // Print the CSRMatrix
    // printCSRMatrix(csrMatrix);

    // Create a PNG image representing the sparsity pattern of the matrix
    createSparsePatternImage(csrMatrix->csr_data, csrMatrix->col_ind, csrMatrix->row_ptr, csrMatrix->num_rows, csrMatrix->num_cols, argv[1]);

    // Assume b is [1, 1, ..., 1]
    double* b = (double*)malloc(csrMatrix->num_rows * sizeof(double));
    double* x = (double*)malloc(csrMatrix->num_rows * sizeof(double)); // Solution vector initialization
    if (b == NULL || x == NULL) {
        fprintf(stderr, "Error allocating memory for b or x vector.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < csrMatrix->num_rows; i++) {
        b[i] = 1.0;
        x[i] = 0.0;
    }

    // Solve via BiCGSTAB and time it
    clock_t start = clock();
    bicgstab(csrMatrix, b, x, 1e-7, 10000);
    
    // Compute BiCGSTAB residual
    double* Ax = (double*)malloc(csrMatrix->num_rows * sizeof(double));
    spmv_csr(csrMatrix, x, Ax);
    double bicgstab_residual = 0.0;
    for (int i = 0; i < csrMatrix->num_rows; i++) {
        bicgstab_residual += (Ax[i] - b[i]) * (Ax[i] - b[i]);
        // printf("%lf ", x[i]); // Use this to print x
    }
    bicgstab_residual = sqrt(bicgstab_residual);

    // Solve via Conjugate Gradient
    conjugate_gradient(csrMatrix, b, x, 1e-7, 10000);
    Ax = (double*)realloc(Ax, csrMatrix->num_rows * sizeof(double));
    spmv_csr(csrMatrix, x, Ax);
    double conj_grad_residual = 0.0;
    for (int i = 0; i < csrMatrix->num_rows; i++) {
        conj_grad_residual += (Ax[i] - b[i]) * (Ax[i] - b[i]);
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