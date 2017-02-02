#define USE_NONLIN_KERNEL 1
#if (USE_NONLIN_KERNEL==1)
void nonlin_kernel(double *vexy,double *dgdxy, double *nonlin,long int length) {
  long int i;

  #pragma simd
  for (i=0;i<length;i++) {
    nonlin[i]=-vexy[i]*dgdxy[length+i]+vexy[length+i]*dgdxy[i];
  }
}
#elif (USE_NONLIN_KERNEL==2)

#include <immintrin.h>
void nonlin_kernel(double *vexy,double *dgdxy, double *nonlin,long int length) {
  long int i;

  __m256d v0,v1,g0,g1;
  __m256d first,second;
  __m256d res;

  for (i=0;i<length;i+=4) {
    v0 = _mm256_loadu_pd(&(vexy[i]));
    v1 = _mm256_loadu_pd(&(vexy[length+i]));
    g0 = _mm256_loadu_pd(&(dgdxy[i]));
    g1 = _mm256_loadu_pd(&(dgdxy[length+i]));

    first  = _mm256_mul_pd(v0,g1);
    second = _mm256_mul_pd(v1,g0);
    res    = _mm256_sub_pd(second,first);
    /*nonlin[i]=-vexy[i]*dgdxy[length+i]+vexy[length+i]*dgdxy[i];*/
    _mm256_storeu_pd(&(nonlin[i]), res);
  }
}
#elif (USE_NONLIN_KERNEL==3)
#include <immintrin.h>
void nonlin_kernel(double *vexy,double *dgdxy, double *nonlin,long int length) {
  long int i;

  __m128d v0,v1,g0,g1;
  __m128d first,second;
  __m128d res;

  for (i=0;i<length;i+=2) {
    v0 = _mm_loadu_pd(&(vexy[i]));
    v1 = _mm_loadu_pd(&(vexy[length+i]));
    g0 = _mm_loadu_pd(&(dgdxy[i]));
    g1 = _mm_loadu_pd(&(dgdxy[length+i]));

    first  = _mm_mul_pd(v0,g1);
    second = _mm_mul_pd(v1,g0);
    res    = _mm_sub_pd(second,first);
    /*nonlin[i]=-vexy[i]*dgdxy[length+i]+vexy[length+i]*dgdxy[i];*/
    _mm_storeu_pd(&(nonlin[i]), res);
  }
}

#endif

/* ======================== zavdxpy ============================== */
/* y(i) = a(i)*x(i)+y(i), with x,y double complex, a double real   */
/* we treat the complex numbers as two real numbers for better
   vectorization. */
#define USE_ZAVDXPY 1

#if (USE_ZAVDXPY==1)
void zavdxpy(const int n, const double *a, const double *x, double *y) {
  int i;

  for (i=0;i<2*n;i++) {
    y[i] = y[i] + a[i/2]*x[i];
  }
}
#elif (USE_ZAVDXPY==2)
void zavdxpy(const int n, const double *a, const double *x, double *y) {
  int i;

  __m256d xv,yv;
  __m256d av,temp;
  __m128d xmm1,xmm2;

  for (i=0;i<2*n;i+=4) {
    /* load x and y into a vector register */
    xv = _mm256_loadu_pd(&(x[i]));
    yv = _mm256_loadu_pd(&(y[i]));

    xmm1 = _mm_loaddup_pd(&(a[i/2]));
    xmm2 = _mm_loaddup_pd(&(a[i/2+1]));
    av =  _mm256_insertf128_pd(av,xmm2, 1);
    av =  _mm256_insertf128_pd(av,xmm1, 0);

    temp = _mm256_mul_pd(av,xv);
    yv = _mm256_add_pd(temp,yv);
    _mm256_storeu_pd(&(y[i]),yv);

    //first  = _mm256_mul_pd(v0,g1);
    //second = _mm256_mul_pd(v1,g0);
    //res    = _mm256_sub_pd(second,first);
    /*nonlin[i]=-vexy[i]*dgdxy[length+i]+vexy[length+i]*dgdxy[i];*/
    //_mm256_storeu_pd(&(nonlin[i]), res);
  }
}
#endif
