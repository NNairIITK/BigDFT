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
#include <iostream>
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



inline
#if DEBUGSPIN
int pthread_spin_lock (const char * str1, int str2, pthread_spinlock_t* spin)
#else
int pthread_spin_lock (const char * , int , pthread_spinlock_t* spin)
#endif
{
#if DEBUGSPIN
  long long int hw1;
  long long int hw2;

  rdtscll(hw1);
  //std::cout<<"locking ("<<pthread_self()<<") : "<<str1<<" "<<str2;

  int err = pthread_spin_lock(spin);

  //std::cout<<"out"<<std::endl;
  rdtscll(hw2);
  if (hw2 - hw1 > 100000)
    // std::cout<<str1<<" "<<str2<<"("<<pthread_self()<<")"<<" : "<<hw2-hw1<<std::endl;
    fprintf (stderr, "%s %d (%d) : %ld\n", str1, str2, (int)pthread_self(),(long int)(hw2-hw1));
  return err;
#else
  return pthread_spin_lock(spin);
#endif
}


inline
#if DEBUGSPIN
int pthread_spin_unlock (const char * str1, int str2, pthread_spinlock_t* spin)
#else
int pthread_spin_unlock (const char * , int , pthread_spinlock_t* spin)
#endif
{
#if DEBUGSPIN
  long long int hw1;
  long long int hw2;

  rdtscll(hw1);
  //std::cout<<"locking ("<<pthread_self()<<") : "<<str1<<" "<<str2;

  int err = pthread_spin_unlock(spin);

  //std::cout<<"out"<<std::endl;
  rdtscll(hw2);
if (hw2 - hw1 > 500)
//  std::cout<<str1<<" "<<str2<<"("<<pthread_self()<<")"<<" : "<<hw2-hw1<<std::endl;
 printf ("%s %d (%d) : %d\n", str1, str2, (int)pthread_self(),(int)(hw2-hw1));
  return err;
#else
  return pthread_spin_unlock(spin);
#endif            
} 
       

template <typename T>
T MAX ( T a, T b)
{
  return ((a < b)? b : a);
}

template <typename T>
T MIN (T a, T b)
{
  return ((a < b)? a : b);
}

template <typename T>
T* alloc_T_random (int size)
{
  T* t = new T[size];
  for (int i=0; i<size; i++)
    t[i] = (T)rand();

  return t;
}

inline int mylog(int x)
{
  int i=0;
  while (x != 0)
    {
      i++;
      x = x>>1;
    }
  return i;
}

inline long long int getcyclecount()
{
  long long int hw1;
  rdtscll(hw1);
  return hw1;
}

#endif
