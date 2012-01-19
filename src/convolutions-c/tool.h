/*  Copyright (C) 2007 Erik Saule, Brice Videau

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
#ifndef PASTEL_TOOL_H
#define PASTEL_TOOL_H
#include <stdlib.h>
typedef    unsigned int      u32;
#include <asm/msr.h>

#ifndef rdtscll

#define rdtscll(t) do { \
     unsigned int __a,__d; \
     asm volatile("rdtsc" : "=a" (__a), "=d" (__d)); \
     (t) = ((unsigned long)__a) | (((unsigned long)__d)<<32); \
} while(0)
//#define rdtscll(val) __asm__ __volatile__("rdtsc" : "=A" (val))
#endif

void rdtsc_(long long unsigned int * t);

#define SECURE_DELETE(X) if (X) {delete X; X=NULL;}
#define SECURE_DELETE_ARRAY(X) if (X) {delete[] X; X=NULL;}


#define ASMPAUSE asm("" : : : "memory")

//X: commande a timer. Y nombre de cycle qui provoque l'affichage
#if DEBUGWHILE
#define TIME_IT(X,Y) {\
  long long int hw1;\
  long long int hw2;\
  rdtscll(hw1);\
  X;\
  rdtscll(hw2);\
  if (hw2 - hw1 > Y)\
     fprintf (stderr, "%s:%d (%d) \"%s\": %ld\n", __FILE__, __LINE__, (int)pthread_self(),#X,(long int)(hw2-hw1));\
}
#else
#define TIME_IT(X,Y) X
#endif

#if DEBUGAUTRE
#define TIME_IT_BEG {\
  long long int hw1;\
  long long int hw2;\
  rdtscll(hw1);

#define TIME_IT_END(X,Y) rdtscll(hw2);\
  if (hw2 - hw1 > Y)\
     fprintf (stderr, "%s:%d (%d) \"%s\" : %ld\n", __FILE__, __LINE__, (int)pthread_self(), X, (long int)(hw2-hw1));\
}
#else
#define TIME_IT_BEG
#define TIME_IT_END(X,Y)
#endif


inline long long int getcyclecount()
{
  long long int hw1;
  rdtscll(hw1);
  return hw1;
}

#endif
