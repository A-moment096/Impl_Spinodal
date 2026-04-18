#ifndef MY_COMPLEX_FFTW_1D
#define MY_COMPLEX_FFTW_1D

#include <math.h>
#include <stddef.h>

typedef double my_complex[2];
static inline void cadd(my_complex a, my_complex b, my_complex out) {
	out[0] = a[0] + b[0];
	out[1] = a[1] + b[1];
}

static inline void csub(my_complex a, my_complex b, my_complex out) {
	out[0] = a[0] - b[0];
	out[1] = a[1] - b[1];
}

static inline void cmul(my_complex a, my_complex b, my_complex out) {
	double re = a[0] * b[0] - a[1] * b[1];
	double im = a[0] * b[1] + a[1] * b[0];
	out[0] = re;
	out[1] = im;
}

static inline void c_exp_pure_image(double theta, my_complex out) {
	out[0] = cos(theta);
	out[1] = sin(theta);
}

static inline void cassign(my_complex dest, my_complex src) {
	dest[0] = src[0];
	dest[1] = src[1];
}

void my_fft_1d_v1(my_complex *in, my_complex *out, size_t N, int sign);
void my_fft_forward_1d_v1(my_complex *in, my_complex *out, size_t N);
void my_fft_backward_1d_v1(my_complex *in, my_complex *out, size_t N);

void my_fft_1d_v2(my_complex *in, my_complex *out, size_t N, int sign);
void my_fft_forward_1d_v2(my_complex *in, my_complex *out, size_t N);
void my_fft_backward_1d_v2(my_complex *in, my_complex *out, size_t N);

#endif