#include "timing.h"



double timing()
{
   struct timeval tp;

   gettimeofday(&tp, NULL);
   return (double) (tp.tv_sec + tp.tv_usec/1000000.0);
}

double timing_() {
   return timing();
}


