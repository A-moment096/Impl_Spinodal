#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "C_fft_1d.h"

#define M_PI 3.14159265358979323846 /* pi */

void my_fft_1d_v1(my_complex *in, my_complex *out, size_t N, int sign)
// N is assumed power of 2
{
    if (N == 1) {
        cassign(out[0], in[0]);
        return;
    }
    my_complex *Ek = (my_complex *)malloc(sizeof(my_complex) * N / 2);
    my_complex *Ek_out = (my_complex *)malloc(sizeof(my_complex) * N / 2);
    my_complex *Ok = (my_complex *)malloc(sizeof(my_complex) * N / 2);
    my_complex *Ok_out = (my_complex *)malloc(sizeof(my_complex) * N / 2);
    if (!Ek || !Ek_out || !Ok || !Ok_out) {
        free(Ek);
        free(Ek_out);
        free(Ok);
        free(Ok_out);
        fprintf(stderr, "%s", "Allocation Failed!");
        return;
    }

    for (size_t i = 0; i < N / 2; ++i) {
        cassign(Ek[i], in[2 * i]);
        cassign(Ok[i], in[2 * i + 1]);
    }
    my_fft_1d_v1(Ek, Ek_out, N / 2, sign);
    my_fft_1d_v1(Ok, Ok_out, N / 2, sign);
    double theta = (double)sign * 2.0 * M_PI / (double)N;
    for (size_t i = 0; i < N / 2; i++) {
        my_complex Wk;
        c_exp_pure_image((double)i * theta, Wk);
        my_complex WOk_k;
        cmul(Wk, Ok_out[i], WOk_k);
        cadd(Ek_out[i], WOk_k, out[i]);
        csub(Ek_out[i], WOk_k, out[i + N / 2]);
    }
    free(Ek);
    free(Ek_out);
    free(Ok);
    free(Ok_out);
    return;
}

void my_fft_1d_v2(my_complex *in, my_complex *out, size_t N, int sign) {
    // reorder
    for (size_t i = 1, j = 0; i < N; i++) {
        size_t bit = N >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            my_complex temp;
            cassign(temp, in[j]);
            cassign(in[j], in[i]);
            cassign(in[i], temp);
        }
    }
    // cassign(out, in);
    memcpy(out, in, sizeof(my_complex) * N);
    for (size_t len = 1; len < N; len <<= 1) {
        my_complex W_len;
        double theta = (double)sign * 2.0 * M_PI / (2.0 * (double)len);
        c_exp_pure_image(theta, W_len);
        for (size_t i = 0; i < N; i += 2 * len) {
            my_complex W = {1.0, 0.0};
            for (size_t j = 0; j < len; j++) {
                size_t k_this = i + j;
                size_t k_next = i + j + len;
                my_complex WO;
                cmul(W, out[k_next], WO);

                my_complex out_this, out_next;
                cadd(out[k_this], WO, out_this);
                csub(out[k_this], WO, out_next);

                cassign(out[k_this], out_this);
                cassign(out[k_next], out_next);
                cmul(W, W_len, W);
            }
        }
    }
}

void my_fft_forward_1d_v1(my_complex *in, my_complex *out, size_t N) {
    my_fft_1d_v1(in, out, N, -1);
}

void my_fft_backward_1d_v1(my_complex *in, my_complex *out, size_t N) {
    my_fft_1d_v1(in, out, N, +1);
    for (size_t i = 0; i < N; i++) {
        out[i][0] /= (double)N;
        out[i][1] /= (double)N;
    }
}

void my_fft_forward_1d_v2(my_complex *in, my_complex *out, size_t N) {
    my_fft_1d_v2(in, out, N, -1);
}

void my_fft_backward_1d_v2(my_complex *in, my_complex *out, size_t N) {
    my_fft_1d_v2(in, out, N, +1);
    for (size_t i = 0; i < N; i++) {
        out[i][0] /= (double)N;
        out[i][1] /= (double)N;
    }
}
