#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <png.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "functions.h"

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix) {
    // Try to open file
    FILE *file = NULL;
    file = fopen(filename, "r"); // Returns pointer to file if successful or NULL if unsuccessful

    // If the file doesn't exist or can't open it, print error message and quit the program
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    int rows = 0, columns = 0, nonzero_values = 0, *row_ind, lower_triangular = 1, upper_triangular = 1;
    double mtx_row, mtx_column, mtx_value;
    char *line = NULL; // Line will point to each line in the file as it's read
    size_t len = 0; // Unsigned int representing the length of the line
    ssize_t read; // Signed int representing the number of characters read. -1 at EOF.

    // Throw away all the lines that start with %, these are comments/headers
    do {
        read = getline(&line, &len, file); // Read a line from the file, do nothing with it
    } while (read != -1 && line[0] == '%'); // Stop when reach EOF or a line that doesn't start with %

    // If reach EOF before finding non-% line
    if (read == -1) {
        fprintf(stderr, "This MTX format is invalid.\n");
        exit(EXIT_FAILURE);
    }

    sscanf(line, "%d %d %d", &rows, &columns, &nonzero_values); // First line of mtx has rows, columns, nonzero values
    free(line);

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

    // Read the values into the CSRMatrix
    int max_col_seen = 0;
    for (int i = 0; i < matrix->num_non_zeros; i++) {
        fscanf(file, "%lf %lf %lf", &mtx_row, &mtx_column, &mtx_value);

        // Assures MTX file is sorted by columns (major), then rows (minor) ascending
        if (mtx_column < max_col_seen) {
            printf("This is improper MTX format. MTX format requires the column-index column to be ascending.\n");
            exit(EXIT_FAILURE);
        } else {
            max_col_seen = mtx_column;
        }

        // Even if 0 values are in the MTX file, they are not stored in the CSRMatrix
        if (mtx_value != 0.0) {
            matrix->csr_data[i] = mtx_value;
            matrix->col_ind[i] = mtx_column - 1; // (1, 1) in MTX is (0, 0) in CSR
            row_ind[i] = mtx_row - 1; // (1, 1) in MTX is (0, 0) in CSR

            // Check if matrix is triangular
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

    // Cumulate row_ptr such that the value at each index is the sum of the values before it
    for (int i = 0; i < rows; i++) {
        matrix->row_ptr[i + 1] += matrix->row_ptr[i];
    }

    // Assign data to MTXRow struct to sort
    // We need to sort because if you don't then CSR will be in the order of the MTX file
    // MTX reads up to down then left to right. CSR needs to read left to right then up to down.
    MTXRow *data = (MTXRow*)malloc(matrix->num_non_zeros * sizeof(MTXRow));
    if (data == NULL) {
        fprintf(stderr, "Error allocating memory for data.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < matrix->num_non_zeros; i++) {
        data[i].row = row_ind[i];
        data[i].col = matrix->col_ind[i];
        data[i].value = matrix->csr_data[i];
    }
    qsort(data, matrix->num_non_zeros, sizeof(MTXRow), compareMTXData); // Sort the rows such to make it read left to right then up to down (Sort by row, then column)

    // Rewrite the sorted data so that it is in the correct order
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

    // Symmetrise the matrix iff user wants to and it is triangular
    if (upper_triangular == 1 || lower_triangular == 1) {
        char symmetrise;
        scanf(" %c", &symmetrise);
        while (symmetrise != 'y' && symmetrise != 'Y' && symmetrise != '1' && symmetrise != 'n' && symmetrise != 'N' && symmetrise != '0') {
            printf("Invalid response. Please enter Y or N: ");
            scanf(" %c", &symmetrise);
        }
        if (symmetrise == 'y' || symmetrise == 'Y' || symmetrise == '1') {
            int og_num_non_zeros = matrix->num_non_zeros;
            // Loop through each non-zero element. If it is not on the diagonal, then add its inverse to the CSR matrix. Inverse being (col, row, value) instead of (row, col, value)
            for (int i = 0; i < og_num_non_zeros; i++) {
                if (matrix->col_ind[i] != row_ind[i]) {
                    matrix->num_non_zeros++; // Adding a new element so increase the number of non-zero elements
                    matrix->col_ind = realloc(matrix->col_ind, matrix->num_non_zeros * sizeof(int));
                    row_ind = realloc(row_ind, matrix->num_non_zeros * sizeof(int));
                    matrix->csr_data = realloc(matrix->csr_data, matrix->num_non_zeros * sizeof(double));
                    matrix->col_ind[matrix->num_non_zeros - 1] = row_ind[i]; // Add the element to the end of the array
                    row_ind[matrix->num_non_zeros - 1] = matrix->col_ind[i]; // Add the element to the end of the array
                    matrix->csr_data[matrix->num_non_zeros - 1] = matrix->csr_data[i]; // Add the element to the end of the array
                }
            }

            // Since all the new elements were added to the end of the array, we need to sort the array again
            MTXRow* data = (MTXRow*)malloc(matrix->num_non_zeros * sizeof(MTXRow));
            if (data == NULL) {
                fprintf(stderr, "Error allocating memory for data.\n");
                exit(EXIT_FAILURE);
            }
            // Assign data to MTXRow struct to sort via qsort
            for (int i = 0; i < matrix->num_non_zeros; i++) {
                data[i].row = row_ind[i];
                data[i].col = matrix->col_ind[i];
                data[i].value = matrix->csr_data[i];
            }
            qsort(data, matrix->num_non_zeros, sizeof(MTXRow), compareMTXData);
            // Rewrite sorted data to existing CSR matrix
            for (int i = 0; i < matrix->num_non_zeros; i++) {
                row_ind[i] = data[i].row;
                matrix->col_ind[i] = data[i].col;
                matrix->csr_data[i] = data[i].value;
            }
            free(data);
            // Re-compute row_ptr
            // Initialize row_ptr to 0
            for (int i = 0; i < rows+1; i++) {
                matrix->row_ptr[i] = 0;
            }
            // Count the number of non-zero elements in each row
            for (int i = 0; i < matrix->num_non_zeros; i++) {
                matrix->row_ptr[row_ind[i] + 1]++;
            }
            // Cumulate row_ptr such that the value at each index is the sum of the values before it
            for (int i = 0; i < rows; i++) {
                matrix->row_ptr[i + 1] += matrix->row_ptr[i];
            }
        }
    }

    free(row_ind);
}

// Sorting criteria for qsort
// CSR format needs to read left to right then up to down, so sort by row then column
int compareMTXData(const void *a, const void *b) {
    // Cast void pointers to MTXRow pointers to access row and col
    MTXRow *rowA = (MTXRow*)a;
    MTXRow *rowB = (MTXRow*)b;

    // Compare rows
    if (rowA->row != rowB->row) {
        return rowA->row - rowB->row;
    } else {
        return rowA->col - rowB->col; // Compare columns if rows are equal
    }
}

// Function to create a PNG image representing the sparsity pattern of a CSR matrix
// This function was adapted from ChatGPT
void createSparsePatternImage(const int *columns, const int *row_ptr, int num_rows, int num_cols, const char *mtx_file_name) {

    // Large matrices would make the image file size too large to view
    if (num_rows > 50000 || num_cols > 50000) {
        fprintf(stderr, "Sorry, the matrix is too large to create a sparsity pattern image.\n");
        return;
    }

    // Remove extension of input matrix file name
    char *dot = strchr(mtx_file_name, '.'); // Points to the . in the file name
    size_t base_name_length = dot != NULL ? dot - mtx_file_name : strlen(mtx_file_name); // Find length of string before the .
    char base_name[256]; // Max file name size = 256
    strncpy(base_name, mtx_file_name, base_name_length); // Copy the string up to the . to base_name
    base_name[base_name_length] = '\0';  // Add EOL char to end of string
    
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
    uint8_t* pixel_data = (uint8_t*)calloc(image_width * image_height, sizeof(uint8_t)); // Intilized all to 0, all pixels are black
    if (pixel_data == NULL) {
        fprintf(stderr, "Error allocating memory for pixel data.\n");
        exit(EXIT_FAILURE);
    }
    // Populate pixel_data based on the sparse matrix
    for (int i = 0; i < num_rows; ++i) {
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; ++j) {
            int col = columns[j];
            // Set the pixel_data to white for each non-zero element
            // Need double for loop to set multiple pixels for each non-zero element due to magnitifcation
            for (int di = 0; di < scale_factor; ++di) {
                for (int dj = 0; dj < scale_factor; ++dj) {
                    pixel_data[(i * scale_factor + di) * image_width + col * scale_factor + dj] = 255; // 255 is white
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

    // Check if the directory exists. If not, create it.
    struct stat st = {0};
    if (stat("Sparsity Pattern Images", &st) == -1) {
        // Create the directory
        #ifdef _WIN32
        _mkdir("Sparsity Pattern Images");
        #else
        mkdir("Sparsity Pattern Images", 0700);
        #endif
    }

    // Write the image data to a file
    char file_name[256];  // Adjust the size as needed
    sprintf(file_name, "Sparsity Pattern Images/%s_sparsity_pattern.png", base_name); // Create file name and assign to file_name
    FILE* fp = fopen(file_name, "wb"); // Open file for writing
    if (!fp) abort(); // If file can't be opened, abort
    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL); // Create png struct to write PNG data
    if (!png) abort(); // If png struct can't be created, abort
    png_infop info = png_create_info_struct(png); // Create png info struct to store image metadata
    if (!info) abort(); // If png info struct can't be created, abort

    if (setjmp(png_jmpbuf(png))) abort(); // If any of the following libpng operations fail, abort 
    png_init_io(png, fp); // Initialize the IO

    // Write image header data
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
    png_write_info(png, info); // Write image info
    png_write_image(png, rows); // Write image data
    png_write_end(png, NULL); // Write the end of the PNG file
    if (png && info)
        png_destroy_write_struct(&png, &info); // Free allocated memory
    if (fp)
        fclose(fp); // Close file

    free(rows);
    free(pixel_data);

    printf("Sparsity pattern image saved to your current directory: %s\n", file_name);
}

// Function to multiply matrix by vector in CSR format. A * x = y
void spmv_csr(const CSRMatrix *A, const double *x, double *y) {
    for (int i = 0; i < A->num_rows; i++) {
        y[i] = 0.0;
        for (int j = A->row_ptr[i]; j < A->row_ptr[i+1]; j++) {
            y[i] += A->csr_data[j] * x[A->col_ind[j]]; // Sum vector value and matrix row value products
        }
    }
}

// Function to solve a linear system of equations using the BiCGSTAB method
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

// Function to solve a linear system of equations using the Conjugate Gradient method
void conjugate_gradient(CSRMatrix* A, double* b, double* x, double tolerance, int max_iterations) {
    double* r = (double*)malloc(A->num_rows * sizeof(double)); // Residual vector
    double* p = (double*)malloc(A->num_rows * sizeof(double)); // Search direction vector
    double* Ap = (double*)malloc(A->num_rows * sizeof(double)); // Temporary storage of A * p
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