#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <pthread.h>

#ifndef HEADER_H
#define HEADER_H

class Args
{
    public:
    int argc = 0;
    int n = 0;
    int m = 0;
    int k = 0;
    int l = 0;
    int r = 0;
    int s = 0;
    char* filename = nullptr;
    double norm = 0;

    double r1_norm1 = 0;
    double r1_norm2 = 0;
    double r2_norm1 = 0;
    double r2_norm2 = 0;
    double r1 = 0;
    double r2 = 0;

    double t1 = 0, t2 = 0, t3 = 0, t4 = 0;

    int t = 0;
    int p = 0;

    double* A = nullptr;
    double* B = nullptr;
    double* X = nullptr;
    
    double* block = nullptr;
    double* T = nullptr;
    double* mat = nullptr;
    double* vec = nullptr;
    bool* idk = nullptr;
    
    pthread_t tid = -1;
    pthread_barrier_t* barrier = nullptr;
    
    double sol_full_time = 0;
    double sol_cpu_time = 0;
    double d_full_time = 0;
    double d_cpu_time = 0;
    int sol_status = 0;

    int init(int n_, int m_, int r_, int s_, char* filename_, int t_, int p_, int argc_, double* A_, double* B_, double* X_, pthread_barrier_t* barrier_)
    {
        n = n_;
        m = m_;
        k = n / m;
        l = n % m;
        r = r_;
        s = s_;
        filename = filename_;
        t = t_;
        p = p_;
        argc = argc_;
        A = A_;
        B = B_;
        X = X_;
        barrier = barrier_;

        block = new double[m * m];
        T = new double[m * 2 * m];
        mat = new double[2 * m * m];
        vec = new double[2 * m];
        idk = new bool[m * m];

        if(!block || !T || !mat || !vec || !idk)
        {
            delete []block;
            delete []T;
            delete []mat;
            delete []vec;
            delete []idk;
            printf("Can't allocate memory in thread %d\n", t); 
            return 1;
        }

        return 0;
    }
    
    void print() const
    {
        printf(
            "thread: %d, n: %d, m: %d, "
            "norm: %lf, p: %d, tid: %ld, cpu_time = %lf, full_time = %lf\n\n",
            t, n, m, norm, p, tid, sol_cpu_time, sol_full_time
        );
    }

    ~Args()
    {
        delete []block;
        delete []T;
        delete []mat;
        delete []vec;
        delete []idk;
    }
};


void GetBlock(double* A, double* block, int n, int m, int p, int q);
void SetBlock(double* A, double* block, int n, int m, int p, int q);
void aboba(double* A, double* mat, int n, int m, int i, int j, int p, int q);
void aboba_set(double* A, double* mat, int n, int m, int i, int j, int p, int q);
void aboba1(double* B, double* vec, int n, int m, int i, int j);
void aboba1_set(double* B, double* vec, int n, int m, int i, int j);
void FillInv(double* A, double* T, int n, int m, int i, int j);
void CalcInv(double* T, int n, int m, int p, int q);
int ReadFromFile(const char* filename, double* matrix, int n);
void InitializeMatrix(double* A, int s, int n);
void PrintMatrix(const double* A, int r, int l, int n);
void CalcB(const double* A, double* B, int n);
void MakeDiag(double* A, Args* a, int n, int m, int i1, int j1, int h, int w);
void ApplyToRight(double* A, Args* a, int n, int m, int i1, int j1, int h, int w, int w1);
void ApplyToVector(double* B, Args* a, int n, int m, int i1, int h, int w);
void FullNullDown(double *A, Args* a, int n, int m, int i1, int j1, int i2, int j2, int h, int w, int h1);
void PartiallyNullDown(double *A, Args* a, int n, int m, int i1, int j1, int i2, int j2, int h, int w, int h1);
void TwiceApplytoRight(double* A, Args* a, int n, int m, int i1, int j1, int i2, int j2, int h, int w, int h1, int w1);
void TwicePartiallyApplytoRight(double* A, Args* a, int n, int m, int i1, int j1, int i2, int j2, int h, int w, int h1, int w1);
void TwiceApplyToVector(double* B, Args* a, int n, int m, int i1, int i2, int h, int w, int h1);
void TwicePartiallyApplyToVector(double* B, Args* a, int n, int m, int i1, int i2, int h, int w, int h1);
double get_full_time(void);
double get_cpu_time(void);
void* thread_func(void* arg);

template <typename T>
int special_reduce_sum(int p, T* a1, T* a2, T* a3, T* a4, int n) {
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static T *r1 = nullptr, *r2 = nullptr, *r3 = nullptr, *r4 = nullptr;
    int i;

    if (p <= 1) return 0;

    pthread_mutex_lock(&m);

    if (r1 == nullptr) {
        r1 = a1;
        r2 = a2;
        r3 = a3;
        r4 = a4;
    } else {

        for (i = 0; i < n; i++) {
            r1[i] += a1[i];
            r2[i] += a2[i];
            r3[i] += a3[i];
            r4[i] += a4[i];
        }
    }

    t_in++;

    if (t_in >= p) {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    } else {
        while (t_in < p) {
            pthread_cond_wait(&c_in, &m);
        }
    }

    if (r1 != a1) {
        for (i = 0; i < n; i++) {
            a1[i] = r1[i];
            a2[i] = r2[i];
            a3[i] = r3[i];
            a4[i] = r4[i];
        }
    }

    t_out++;

    if (t_out >= p) {
        t_in = 0;
        r1 = r2 = r3 = r4 = nullptr;
        pthread_cond_broadcast(&c_out);
    } else {
        while (t_out < p) {
            pthread_cond_wait(&c_out, &m);
        }
    }

    pthread_mutex_unlock(&m);
    return 0;
}

template <typename T>
int reduce_sum(int p, T* a, int n) {
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static T *r = nullptr;
    int i;

    if (p <= 1) return 0;

    pthread_mutex_lock(&m);

    if (r == nullptr) {
        r = a;
    } else {

        for (i = 0; i < n; i++) {
            r[i] += a[i];
        }
    }

    t_in++;

    if (t_in >= p) {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    } else {
        while (t_in < p) {
            pthread_cond_wait(&c_in, &m);
        }
    }

    if (r != a) {
        for (i = 0; i < n; i++) {
            a[i] = r[i];

        }
    }

    t_out++;

    if (t_out >= p) {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    } else {
        while (t_out < p) {
            pthread_cond_wait(&c_out, &m);
        }
    }

    pthread_mutex_unlock(&m);
    return 0;
}

#endif 
    
