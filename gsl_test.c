/*
gsl_test

This program tests your installation of gsl and libgsl-dev to make sure you can use
the excellent gsl random-number generators.  Here we demonstrate using the 
"Mersenne Twister" pseudorandom number generator.

compile this via

gcc -o gsl_test gsl_test.c -lgsl

Cameron Abrams -- cfa22@drexel.edu
*/
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main ( int argc, char * argv[] ) {
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,5488);
  int i;
  float sigma = 0.25;
  float mu = 4.0;

  printf("Here are some random floats uniformly distributed on [0,1]:\n");
  for (i=0;i<10;i++) {
    printf("%.3f\n",gsl_rng_uniform(r));
  }

  printf("Here are some random ints uniformly picked from [0,%d]:\n",10);
  for (i=0;i<10;i++) {
    printf("%ld\n",gsl_rng_uniform_int(r,10));
  }

  printf("Here are some floats drawn randomly from a Gaussian distribution with sigma=%.5f:\n",sigma);
  for (i=0;i<10;i++) {
    printf("%.3f\n",gsl_ran_gaussian(r,sigma));
  }

  printf("Here are a few integers drawn randomly from a Poisson distribution with mean=%.5f:\n",mu);
  for (i=0;i<10;i++) {
    printf("%d\n",gsl_ran_poisson(r,mu));
  }
  printf("That's all!\n");
  gsl_rng_free(r);
}
