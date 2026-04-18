#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "C_fft_1d.h"

int main() {
#define N 8
    my_complex test[N];
    for (double i = 0.0; i < (double)N; i += 1.0) {
        test[(int)i][0] = i;
        test[(int)i][1] = 0.0;
    }
    my_complex *out_v1 = (my_complex *)malloc(sizeof(my_complex) * N);
    my_complex *out_v2 = (my_complex *)malloc(sizeof(my_complex) * N);

    if (!out_v1 || !out_v2) {
        free(out_v1);
        free(out_v2);
    }

    my_fft_1d_v1(test, out_v1, N, -1);
    my_fft_1d_v2(test, out_v2, N, -1);

    for (int i = 0; i < N; i++) {
        printf("%f, %f,\t|\t%f, %f\n", out_v1[i][0], out_v1[i][1], out_v2[i][0], out_v2[i][1]);
    }

    printf("%s\n", "------------------------------------");

    my_fft_backward_1d_v1(out_v1, out_v1, N);
    my_fft_backward_1d_v2(out_v2, out_v2, N);

    for (int i = 0; i < N; i++) {
        printf("%f, %f,\t|\t%f, %f\n", out_v1[i][0], out_v1[i][1], out_v2[i][0], out_v2[i][1]);
    }

    free(out_v1);
    free(out_v2);
}