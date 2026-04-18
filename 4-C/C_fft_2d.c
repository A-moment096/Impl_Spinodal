#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "C_fft_1d.h"

typedef double my_complex[2];

// Transform between row-major to column-major
void my_fft_util_transpose(
    my_complex *in,
    my_complex *out,
    size_t N0,
    size_t N1) {

    size_t N_full = N0 * N1;
    my_complex *transformed;
    char inplace = (in == out);
    if (inplace) {
        transformed = (my_complex *)malloc(sizeof(my_complex) * N_full);
        if (!transformed) {
            free(transformed);
            return;
        }
    } else {
        transformed = out;
    }

    for (size_t j = 0; j < N1; j++) {
        for (size_t i = 0; i < N0; i++) {
            cassign(transformed[j + i * N1], in[i + j * N0]);
        }
    }

    if (inplace) {
        memcpy(out, transformed, sizeof(my_complex) * N_full);
        free(transformed);
        return;
    }
    return;
}

void my_fft_2d_v1(
    my_complex *in,
    my_complex *out,
    size_t N0,
    size_t N1,
    int sign) {
    // N0 and N1 must be power of 2
    // row major: [i,j] = j+i*N0
    size_t N_full = N0 * N1;

    // transform row first
    my_complex *row_transformed_in = (my_complex *)malloc(sizeof(my_complex) * N_full);
    if (!row_transformed_in) {
        free(row_transformed_in);
        return;
    }
    for (size_t i = 0; i < N_full; i += N0) {
        my_fft_1d_v1(&(in[i]), &(row_transformed_in[i]), N0, -1);
    }
    // then transform column

    // perform inplace transpose from row-major to column-major
    my_fft_util_transpose(row_transformed_in, row_transformed_in, N0, N1);

    for (size_t i = 0; i < N_full; i += N1) {
        my_fft_1d_v1(&(row_transformed_in[i]), &(out[i]), N1, -1);
    }

    my_fft_util_transpose(out, out, N1, N0);

    // temp input column
    // my_complex *column_in = (my_complex *)malloc(sizeof(my_complex) * N1);
    // // temp output column
    // my_complex *column_out = (my_complex *)malloc(sizeof(my_complex) * N1);
    // if (!column_out || !column_in) {
    //     free(column_in);
    //     free(column_out);
    //     return;
    // }

    free(row_transformed_in);
    // free(column_in);
    // free(column_out);
    return;
}

int test_my_fft_util_transpose() {
    /*
    0 1
    2 3
    4 5
    |
    V
    0 2 4
    1 3 5
    */
    size_t N0 = 2, N1 = 3;
    size_t N = N0 * N1;
    my_complex *in_trans = (my_complex *)malloc(sizeof(my_complex) * N);
    if (!in_trans) {
        free(in_trans);
        return -1;
    }
    for (int i = 0; i < N; i++) {
        in_trans[i][0] = (double)i;
        in_trans[i][1] = 0.0;
    }
    for (int j = 0; j < N1; j++) {
        for (int i = 0; i < N0; i++) {
            printf("(%.3f,%.3f) ", in_trans[i + j * N0][0], in_trans[i + j * N0][1]);
        }
        printf("\n");
    }
    printf("\n");
    my_fft_util_transpose(in_trans, in_trans, N0, N1);
    for (int j = 0; j < N0; j++) {
        for (int i = 0; i < N1; i++) {
            printf("(%.3f,%.3f) ", in_trans[i + j * N1][0], in_trans[i + j * N1][1]);
        }
        printf("\n");
    }
    free(in_trans);
    return 0;
}

int my_print_complex_matrix(my_complex *matrix, size_t N0, size_t N1) {
    for (int j = 0; j < N0; j++) {
        for (int i = 0; i < N1; i++) {
            printf("(%.3f,%.3f) ", matrix[i + j * N0][0], matrix[i + j * N0][1]);
        }
        printf("\n");
    }
}

int test_my_fft_2d_v1() {
#define N0 4
#define N1 4
    size_t N = N0 * N1;
    my_complex *in = (my_complex *)malloc(sizeof(my_complex) * N);
    if (!in) {
        free(in);
        return -1;
    }
    for (int i = 0; i < N; i++) {
        in[i][0] = (double)(1 << i);
        in[i][1] = 0.0;
    }

    my_fft_2d_v1(in, in, N0, N1, -1);
    my_print_complex_matrix(in, N0, N1);

    //
    free(in);
    return 0;
}

int main() {
    return test_my_fft_2d_v1();
}
