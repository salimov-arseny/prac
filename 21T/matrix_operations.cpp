#include "header.h"

void aboba1(double* B, double* vec, int n, int m, int i, int j)
{
    //double tmp = clock();
    int k = n / m;
    int l = n % m;

    int h = (i < k) ? m : l;

    int h1 = (j < k) ? m : l;

    int s;
    
    for(s = 0; s < h1; s++)
        vec[s] = B[j * m + s];

    for(s = 0; s < h; s++)
        vec[s + m] = B[i * m + s];
    //printf("aboba1() %lf\n", (clock() - tmp)/CLOCKS_PER_SEC);
}

void aboba1_set(double* B, double* vec, int n, int m, int i, int j)
{
    //double tmp = clock();
    int k = n / m;
    int l = n % m;

    int h = (i < k) ? m : l;

    int h1 = (j < k) ? m : l;

    int s;
    
    for(s = 0; s < h1; s++)
        B[j * m + s] = vec[s];
    
    for(s = 0; s < h; s++)
        B[i * m + s] = vec[s + m];
    //printf("aboba1() %lf\n", (clock() - tmp)/CLOCKS_PER_SEC);
}

void FillInv(double* A, double* T, int n, int m, int p, int q)
{
    int k = n / m;
    int l = n % m;

    int h = p < k ? m : l;

    int cnt = 0;
    int cnt1 = 0;

    double* start = A + p * (k * m * m + m * l) + q * m;

    for(int i = 0; i < h; i++)
        for(int j = 0; j < h; j++)
        {
            T[i * 2 * m + j] = start[cnt + cnt1];
            T[i * 2 * m + j + m] = (i == j) ? 1 : 0; 
            cnt1 += (cnt + 1) % h == 0 ? n - h : 0;
            cnt++;
        }
}

void CalcInv(double* T, int n, int m, int p, int q)
{
    int k = n / m;
    int l = n % m;
    int h = p < k ? m : l;
    double sum = 0;
    (void)q;

    for(int i = 0; i < h; i++)
        T[i * 2 * m + i + m] = 1 / T[i * 2 * m + i];
    
    for(int i = h - 2; i >= 0; i--)
        for(int j = i + 1; j < h; j++)
        {
            sum = 0;
            for(int k = i + 1; k < j + 1; k++)
                sum += T[i * 2 * m + k] * T[k * 2 * m + j + m];
            
            T[i * 2 * m + j + m] = -sum/T[i * 2 * m + i];
        }
}

int ReadFromFile(const char* filename, double* A, int n) 
{
    std::ifstream f(filename);
    int cnt = 0;
    if (!f) 
    {
        printf("Error: Unable to open file %s\n", filename);
        return 1;
    }

    for (int i = 0; i < n; ++i) 
        for (int j = 0; j < n; ++j, cnt++)        
            if (!(f >> A[i * n + j])) 
            {
                printf("Error: Invalid matrix format in file %s\n", filename);
                f.close();
                return 1;
            }

    if(cnt != n*n)
    {
        printf("Error: Not enough values in file %s\n", filename);
        f.close();
        return 1;
    }

    f.close();
    return 0;
}

void InitializeMatrix(double* A, int s, int n)
{
    auto f = [](int s, int i, int j, int n) -> double
    {
        switch(s)
        {
            case 1: return n - std::max(i, j);
            case 2: return std::max(i, j) + 1;
            case 3: return std::abs(i - j);
            case 4: return 1.0 / (i + j + 1);
            default: return 0.0;
        }
    };

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            A[i * n + j] = f(s, i, j, n);
}

void PrintMatrix(const double* A, int r, int l, int n) 
{
    for (int i = 0; i < std::min(r, l); ++i) 
    {
        for (int j = 0; j < std::min(r, n); ++j) 
        {
            printf(" %10.3e", A[i * n + j]);
            //if((j + 1) % 4 == 0)
              //  printf("  ");
        }

        //if((i + 1) % 4 == 0)
          //  printf("\n");
        printf("\n");
    }
}

void CalcB(const double* A, double* B, int n) 
{
    for (int i = 0; i < n; i++) 
    {
        B[i] = 0;
        for (int k = 0; k <= (n + 1) / 2; k++)  
            if (2 * k < n) 
                B[i] += A[i * n + (2 * k)];
    }
}

void GetBlock(double* A, double* block, int n, int m, int p, int q)
{
    int l = n % m;
    int k = n / m;
    int i, j;
    int cnt = 0;
    int cnt1 = 0;

    //double tmp = clock();

    double *start = A + p *(k * m * m + m * l) + q * m;
    if(p < k)
    {
        if(q < k)
        {
            for(i = 0; i < m; i++)
                for(j = 0; j < m; j++)
                {
                    block[i * m + j] = start[cnt + cnt1];
                    cnt1 += (cnt + 1) % m == 0 ? n-m : 0;
                    cnt++;
                }
        }
        else
        {
            for(i = 0; i < m; i++)
                for(j = 0; j < l; j++)
                {
                    block[i * m + j] = start[cnt + cnt1];
                    cnt1 += (cnt + 1) % l == 0 ? n-l : 0;
                    cnt++;
                }
        }
    }
    else 
    {
        if(q < k)
        {
            for(i = 0; i < l; i++)
                for(j = 0; j < m; j++)
                {
                    block[i * m + j] = start[cnt + cnt1];
                    cnt1 += (cnt + 1) % m == 0 ? n-m : 0;
                    cnt++;
                }
        }
        else
        {
            for(i = 0; i < l; i++)
                for(j = 0; j < l; j++)
                {
                    block[i * m + j] = start[cnt + cnt1];
                    cnt1 += (cnt + 1) % l == 0 ? n-l : 0;
                    cnt++;
                }
        }
    }
    //printf("get_block() %lf\n", (clock() - tmp)/CLOCKS_PER_SEC);
}

void SetBlock(double* A, double* block, int n, int m, int p, int q)
{
    int l = n % m;
    int k = n / m;
    int i, j;
    int cnt = 0;
    int cnt1 = 0;

    //double tmp = clock();

    double *start = A + p *(k * m * m + m * l) + q * m;
    if(p < k)
    {
        if(q < k)
        {
            for(i = 0; i < m; i++)
                for(j = 0; j < m; j++)
                {
                    start[cnt + cnt1] = block[i * m + j];
                    cnt1 += (cnt + 1) % m == 0 ? n-m : 0;
                    cnt++;
                }
        }
        else
        {
            for(i = 0; i < m; i++)
                for(j = 0; j < l; j++)
                {
                    start[cnt + cnt1] = block[i * m + j];
                    cnt1 += (cnt + 1) % l == 0 ? n-l : 0;
                    cnt++;
                }
        }
    }
    else 
    {
        if(q < k)
        {
            for(i = 0; i < l; i++)
                for(j = 0; j < m; j++)
                {
                    start[cnt + cnt1] = block[i * m + j];
                    cnt1 += (cnt + 1) % m == 0 ? n-m : 0;
                    cnt++;
                }
        }
        else
        {
            for(i = 0; i < l; i++)
                for(j = 0; j < l; j++)
                {
                    start[cnt + cnt1] = block[i * m + j];
                    cnt1 += (cnt + 1) % l == 0 ? n-l : 0;
                    cnt++;
                }
        }
    }
    //printf("set_block() %lf\n", (clock() - tmp)/CLOCKS_PER_SEC);
}

void MakeDiag(double* A, Args* a, int n, int m, int i1, int j1, int h, int w)
{
    double t = get_full_time();
    double c1, c2, cos, sin, sq, eps = 1e-15 * a->norm;
    GetBlock(A, a->block, n, m, i1, j1);

    for(int j = 0; j < w; j++)
    {
        for(int i = j + 1; i < h; i++)
        {
            c1 = a->block[j * m + j];
            c2 = a->block[i * m + j];
            
            sq = sqrt(c1 * c1  + c2 * c2);
            //printf("i = %d j = %d c1 = %10.3e c2 = %10.3e\n", i, j, c1, c2);
            //PrintMatrix(block, h, h, h);
            if(sq > eps)
            {
                cos = c1 / sq;
                sin = -c2 / sq;
                a->T[i * 2 * m + j] = cos;
                a->T[i * 2 * m + j + m] = sin;

                a->idk[i * m + j] = true;

                a->block[i * m + j] = 0;
                a->block[j * m + j] = sq;
                for(int q = j + 1; q < w; q++) // j + 1
                {
                    c1 = a->block[i * m + q]; 
                    c2 = a->block[j * m + q]; 
                    a->block[i * m + q] = cos * c1 + sin * c2;
                    a->block[j * m + q] = cos * c2 - sin * c1;
                }
            }
            else
            {
                a->block[i * m + j] = 0;
                a->block[j * m + j] = sq;
                a->idk[i * m + j] = false;
            }
        }
        
    }

    SetBlock(A, a->block, n, m, i1, j1);

    t = get_full_time() - t;
    a->t1 += t;
}

void ApplyToRight(double* A, Args* a, int n, int m, int i1, int j1, int h, int w, int w1)
{
    double t = get_full_time();
    double c1, c2, sin, cos;

    GetBlock(A, a->block, n, m, i1, j1);

    for(int j = 0; j < w; j++)
    {
        for(int i = j + 1; i < h; i++)
        {
            if(a->idk[i * m + j] == true)
            {
                for(int q = 0; q < w1; q++)
                {
                    c1 = a->block[i * m + q];
                    c2 = a->block[j * m + q];
                    cos = a->T[i * 2 * m + j]; 
                    sin = a->T[i * 2 * m + j + m];
                    a->block[i * m + q] = c1 * cos + c2 * sin;
                    a->block[j * m + q] = c2 * cos - c1 * sin;
                }
            }
        }
    }
    SetBlock(A, a->block, n, m, i1, j1);

    t = get_full_time() - t;
    a->t1 += t;
}

void ApplyToVector(double* B, Args* a, int n, int m, int i1, int h, int w)
{
    double t = get_full_time();
    double c1, c2, sin, cos;
    (void)n;

    for(int j = 0; j < w; j++)
    {
        for(int i = j + 1; i < h; i++)
        {
            if(a->idk[i * m + j] == true)
            {
                c1 = B[i + i1 * m];
                c2 = B[j + i1 * m];
                //printf("i = %d j = %d c1 = %10.3e c2 = %10.3e\n", i, j, c1, c2);
                cos = a->T[i * 2 * m + j];
                sin = a->T[i * 2 * m + j + m];
                B[i + i1 * m] = cos * c1 + sin * c2;
                B[j + i1 * m] = cos * c2 - sin * c1;
            }
        }
    }

    t = get_full_time() - t;
    a->t1 += t;
}

void aboba(double* A, double* mat, int n, int m, int i, int j, int p, int q)
{
    int k = n / m;
    int l = n % m;

    int cnt1 = 0;
    int cnt = 0;

    int h = (i < k) ? m : l;
    int w = (j < k) ? m : l;

    int h1 = (p < k) ? m : l;
    int w1 = (q < k) ? m : l;

    int t, r;

    //double tmp = clock();

    double* start = A + p * (k * m * m + l * m) + q * m;

        for(t = 0; t < h1; t++)
            for(r = 0; r < w1; r++)
            {
                mat[t * m + r] = start[cnt + cnt1];
                cnt1 += (cnt + 1) % w1 == 0 ? n-w1 : 0;
                cnt++;
            }
    

    cnt1 = cnt = 0;
    start = A + i * (k * m * m + l * m) + j * m;

    for(t = 0; t < h; t++)
        for(r = 0; r < w; r++)
        {
            mat[t * m + r + m * m] = start[cnt + cnt1];
            cnt1 += (cnt + 1) % w == 0 ? n-w : 0;
            cnt++;
        }
    //printf("aboba() %lf\n", (clock() - tmp)/CLOCKS_PER_SEC);
}

void aboba_set(double* A, double* mat, int n, int m, int i, int j, int p, int q)
{
    int k = n / m;
    int l = n % m;

    int cnt1 = 0;
    int cnt = 0;

    int h = (i < k) ? m : l;
    int w = (j < k) ? m : l;

    int h1 = (p < k) ? m : l;
    int w1 = (q < k) ? m : l;

    int t, r;

    //double tmp = clock();

    double* start = A + p * (k * m * m + l * m) + q * m;

    for(t = 0; t < h1; t++)
        for(r = 0; r < w1; r++)
        {
            start[cnt + cnt1] = mat[t * m + r];
            cnt1 += (cnt + 1) % w1 == 0 ? n-w1 : 0;
            cnt++;
        }

    cnt1 = cnt = 0;
    start = A + i * (k * m * m + l * m) + j * m;

    for(t = 0; t < h; t++)
        for(r = 0; r < w; r++)
        {
            start[cnt + cnt1] = mat[t * m + r + m * m] ;
            cnt1 += (cnt + 1) % w == 0 ? n-w : 0;
            cnt++;
        }
    //printf("aboba_set() %lf\n", (clock() - tmp)/CLOCKS_PER_SEC);
}

void FullNullDown(double *A, Args* a, int n, int m, int i1, int j1, int i2, int j2, int h, int w, int h1)
{
    double t = get_full_time();
    double c1, c2, sin, cos, sq, eps = 1e-15 * a->norm;

    (void)h;

    aboba(A, a->mat, n, m, i1, j1, i2, j2);

    for(int j = 0; j < w; j++)  
    {
        for(int i = 0; i < h1; i++)
        {
            c1 = a->mat[j * m + j];
            c2 = a->mat[i * m + j + m * m];
            //printf("v = %d s = %d i = %d j = %d c1 = %10.3e c2 = %10.3e\n", v, s, i, j, c1, c2);
            //PrintMatrix(mat, h + 1, h + 1, m);
            sq = sqrt(c1 * c1 + c2 * c2);
            if(sq > eps)
            {
                a->idk[i * m + j] = true;
                cos = c1 / sq;
                sin = -c2 / sq;
                a->T[i * 2 * m + j] = cos;
                a->T[i * 2 * m + j + m] = sin;

                a->mat[i * m + j + m * m] = 0;
                a->mat[j * m + j] = sq;

                for(int q = j + 1; q < w; q++) 
                {
                    c1 = a->mat[j * m + q];
                    c2 = a->mat[i * m + q + m * m];    
                    a->mat[j * m + q] = cos * c1 - sin * c2;
                    a->mat[i * m + q + m * m] = cos * c2 + sin * c1;
                }
            }
            else
            {
                a->idk[i * m + j] = false;
                a->mat[i * m + j + m * m] = 0;
                a->mat[j * m + j] = sq;
            }
        }
    }
    aboba_set(A, a->mat, n, m, i1, j1, i2, j2);

    t = get_full_time() - t;
    a->t2 += t;
}

void PartiallyNullDown(double *A, Args* a, int n, int m, int i1, int j1, int i2, int j2, int h, int w, int h1)
{
    double t = get_full_time();

    double c1, c2, sin, cos, sq, eps = 1e-15 * a->norm;

    (void)h;

    aboba(A, a->mat, n, m, i1, j1, i2, j2);

    for(int j = 0; j < w; j++)
    {
        for(int i = 0; i < std::min(j + 1, h1); i++)
        {
            c1 = a->mat[j * m + j];
            c2 = a->mat[i * m + j + m * m];
            //printf("v = %d s = %d i = %d j = %d c1 = %10.3e c2 = %10.3e\n", v, s, i, j, c1, c2);
            //PrintMatrix(mat, h + 1, h + 1, m);
            sq = sqrt(c1 * c1 + c2 * c2);
            if(sq > eps)
            {
                a->idk[i * m + j] = true;
                cos = c1 / sq;
                sin = -c2 / sq;
                a->T[i * 2 * m + j] = cos;
                a->T[i * 2 * m + j + m] = sin;

                a->mat[i * m + j + m * m] = 0;
                a->mat[j * m + j] = sq;

                for(int q = j + 1; q < w; q++) 
                {
                    c1 = a->mat[j * m + q];
                    c2 = a->mat[i * m + q + m * m];    
                    a->mat[j * m + q] = cos * c1 - sin * c2;
                    a->mat[i * m + q + m * m] = cos * c2 + sin * c1;
                }
            }
            else
            {
                a->idk[i * m + j] = false;
                a->mat[i * m + j + m * m] = 0;
                a->mat[j * m + j] = sq;
            }
        }
    }

    aboba_set(A, a->mat, n, m, i1, j1, i2, j2);
    t = get_full_time() - t;
    a->t3 += t;
}

void TwiceApplytoRight(double* A, Args* a, int n, int m, int i1, int j1, int i2, int j2, int h, int w, int h1, int w1)
{
    double t = get_full_time();
    double c1, c2, sin, cos;

    (void)h;

    aboba(A, a->mat, n, m, i1, j1, i2, j2);

    for(int j = 0; j < w; j++)
    {
        for(int i = 0; i < h1; i++)
        {
            if(a->idk[i * m + j] == true)
            {
                cos = a->T[i * 2 * m + j];
                sin = a->T[i * 2 * m + j + m];
                for(int q = 0; q < w1; q++)
                {
                    c1 = a->mat[j * m + q];
                    c2 = a->mat[i * m + q + m * m];
                    a->mat[j * m + q] = cos * c1 - sin * c2;
                    a->mat[i * m + q + m * m] = cos * c2 + sin * c1; 
                }
            }
        }
    }
    aboba_set(A, a->mat, n, m, i1, j1, i2, j2);
    t = get_full_time() - t;

    a->t2 += t;
}

void TwicePartiallyApplytoRight(double* A, Args* a, int n, int m, int i1, int j1, int i2, int j2, int h, int w, int h1, int w1)
{
    double t = get_full_time();
    double c1, c2, sin, cos;

    (void)h;

    aboba(A, a->mat, n, m, i1, j1, i2, j2);

    for(int j = 0; j < w; j++)
    {
        for(int i = 0; i < std::min(j + 1, h1); i++)
        {
            if(a->idk[i * m + j] == true)
            {
                cos = a->T[i * 2 * m + j];
                sin = a->T[i * 2 * m + j + m];
                for(int q = 0; q < w1; q++)
                {
                    c1 = a->mat[j * m + q];
                    c2 = a->mat[i * m + q + m * m];
                    a->mat[j * m + q] = cos * c1 - sin * c2;
                    a->mat[i * m + q + m * m] = cos * c2 + sin * c1; 
                }
            }
        }
    }
    aboba_set(A, a->mat, n, m, i1, j1, i2, j2);
    t = get_full_time() - t;
    a->t3 += t;
}

void TwiceApplyToVector(double* B, Args* a, int n, int m, int i1, int i2, int h, int w, int h1)
{
    double t = get_full_time();
    double c1, c2, sin, cos;

    (void)h;

    aboba1(B, a->vec, n, m, i1, i2);

    for(int j = 0; j < w; j++)
    {
        for(int i = 0; i < h1; i++)
        {
            if(a->idk[i * m + j] == true)
            {
                c1 = a->vec[j];
                c2 = a->vec[i + m];

                cos = a->T[i * 2 * m + j];
                sin = a->T[i * 2 * m + j + m];

                a->vec[j] = cos * c1 - sin * c2;
                a->vec[i + m] = cos * c2 + sin * c1;
            }
        }
    }

    aboba1_set(B, a->vec, n, m, i1, i2);

    t = get_full_time() - t;
    a->t2 += t;
}

void TwicePartiallyApplyToVector(double* B, Args* a, int n, int m, int i1, int i2, int h, int w, int h1)
{
    double t = get_full_time();
    double c1, c2, sin, cos;

    (void)h;

    aboba1(B, a->vec, n, m, i1, i2);

    for(int j = 0; j < w; j++)
    {
        for(int i = 0; i < std::min(j + 1, h1); i++)
        {
            if(a->idk[i * m + j] == true)
            {
                c1 = a->vec[j];
                c2 = a->vec[i + m];

                cos = a->T[i * 2 * m + j];
                sin = a->T[i * 2 * m + j + m];

                a->vec[j] = cos * c1 - sin * c2;
                a->vec[i + m] = cos * c2 + sin * c1;
            }
        }
    }

    aboba1_set(B, a->vec, n, m, i1, i2);
    t = get_full_time() - t;
    a->t3 += t;
}

double get_full_time() 
{
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec/1e6;
}

double get_cpu_time() 
{
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec/1e6;
}

void* thread_func(void* arg)
{
    Args* a = (Args*)arg;

    double* A = a->A;
    double* B = a->B;
    double* X = a->X;
    
    int nproc = get_nprocs();
    //printf("%d\n", nproc);
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    CPU_SET(nproc - 1 - (a->t % (nproc)), &cpu);
    pthread_setaffinity_np(pthread_self(), sizeof(cpu), &cpu);

    int n = a->n;
    int m = a->m;

    int k = n / m; 
    int l = n % m; 

    int t = (l == 0) ? k : k + 1; //число шагов

    int h, w, h1, w1, s, i, j, q, v, z, x, y, hp, dg, dob, g, f, r, a2, b2, pwr, i1;

    double c1 = 0;

    dg = 2;
    dob = 1;

    a2 = (int)(log(a->p)/log(2));
    pwr = int(pow(2, a2));
    b2 = a->p - pwr;
    
    r = t / a->p; 
    f = t % a->p; //высота "остаточной" матрицы
    
    hp = (a->t < f) ? r + 1 : r; //высота "своей" матрицы 

    for(q = 0; q < hp; q++)
    {
        h = (q * a->p + a->t < k) ? m : l;
        memset(A + (q * a->p + a->t) * n * m, 0, n * h * sizeof(double));
        memset(B + (q * a->p + a->t) * m, 0, h * sizeof(double));
    }

    pthread_barrier_wait(a->barrier);

    memset(a->block, 0, m * m * sizeof(double));
    memset(a->T, 0, m * 2 * m * sizeof(double));
    memset(a->mat, 0, m * 2 * m * sizeof(double));
    memset(a->vec, 0, m * 2 * sizeof(double));
    memset(a->idk, false, m * m * sizeof(bool));

    if(a->t == 0)
    {
        if (a->s == 0 && a->argc == 7) 
        {
            if(ReadFromFile(a->filename, A, n) != 0)
            {
                a->sol_status = -1;
            }
        }
        else
            InitializeMatrix(A, a->s, n);
    }

    reduce_sum(a->p, &a->sol_status, 1);

    if(a->sol_status < 0)
    {
        return nullptr;
    }

    if(a->t == 0)
    {
        CalcB(A, B, a->n);

        a->norm = std::abs(A[0]);
        for(i = 0; i < n; i++)
        {
            for(j = 0; j < n; j++)
            {
                c1 = std::abs(A[i * n + j]);
                if(c1 < a->norm)
                    a->norm = c1;
            }
        }

        printf("Matrix A:\n");

        PrintMatrix(A, a->r, a->n, a->n);
    }

    reduce_sum(a->p, &a->norm, 1);

    a->sol_cpu_time = get_cpu_time();
    a->sol_full_time = get_full_time();

    for(s = 0; s < t; s++)
    {
        x = s / a->p;
        y = s % a->p;
        g = (a->t < y) ? x + 1 : x; //g * a->p + a->t - главный блок у потока a->t 

        h = (g * a->p + a->t < k) ? m : l;
        w = (s < k) ? m : l;

        if(h == w && g * a->p + a->t < t)
        {
            //printf("Приводим (%d, %d) к диагональному\n", g * a->p + a->t, s);
            MakeDiag(A, a, n, m, g * a->p + a->t, s, h, w);

            for(v = s + 1; v < t; v++)
            {
                w1 = (v < k) ? m : l;
                //printf("Применяем к (%d, %d)\n", g * a->p + a->t, v);
                ApplyToRight(A, a, n, m, g * a->p + a->t, v, h, w, w1);
            }

            ApplyToVector(B, a, n, m, g * a->p + a->t, h, w);
            
            //printf("ГООООООООООООООООООООЛ %d at step %d\n", a->t, s);

            for(v = g + 1; v < hp; v++)
            {

                //printf("Зануляем (%d, %d) \n", v * a->p + a->t, s);
                h1 = (v*a->p + a->t < k) ? m : l;

                FullNullDown(A, a, n, m, v * a->p + a->t, s, g * a->p + a->t, s, h, w, h1);
                
                for(z = s + 1; z < t; z++)
                {
                    w1 = z < k ? m : l;
                    //printf("Применяем к (%d, %d) \n", v * a->p + a->t, z);
                    TwiceApplytoRight(A, a, n, m, v*a->p + a->t, z, g*a->p + a->t, z, h, w, h1, w1);
                }

                TwiceApplyToVector(B, a, n, m, v*a->p + a->t, g*a->p + a->t, h, w, h1);
            }
        }

        pthread_barrier_wait(a->barrier);
        /*if(a->t == 0)
        {
            printf("step %d\n\n\n\n", s);
            
            PrintMatrix(A, a->r, n, n);
            printf("\n\n\n\n");
        }

        pthread_barrier_wait(a->barrier);*/

        if((a->t >= y && a->t < y + b2 && pwr + g * a->p + a->t < t) || (a->t < y && a->t + a->p < y + b2 && pwr + g * a->p + a->t < t))
        {
            //printf("Зануляем (%d, %d)\n", pwr + g * a->p + a->t, s);
            //aboba(A, a->mat, n, m, pwr + g * a->p + a->t, s, g * a->p + a->t, s);
            h = (g * a->p + a->t < k) ? m : l;
            w = (s < k) ? m : l;
            h1 = (pwr + g * a->p + a->t < k) ? m : l;

            if(h1 == w)
            {
                PartiallyNullDown(A, a, n, m, pwr + g * a->p + a->t, s, g * a->p + a->t, s, h, w, h1);

                for(z = s + 1; z < t; z++)
                {
                    w1 = z < k ? m : l;
                    TwicePartiallyApplytoRight(A, a, n, m, pwr + g * a->p + a->t, z, g * a->p + a->t, z, h, w, h1, w1);
                }

                TwicePartiallyApplyToVector(B, a, n, m, pwr + g * a->p + a->t, g * a->p + a->t, h, w, h1);
            }
            else
            {   
                FullNullDown(A, a, n, m, pwr + g * a->p + a->t, s, g * a->p + a->t, s, h, w, h1);

                for(z = s + 1; z < t; z++)
                {
                    w1 = z < k ? m : l;
                    TwiceApplytoRight(A, a, n, m, pwr + g * a->p + a->t, z, g * a->p + a->t, z, h, w, h1, w1);
                }

                TwiceApplyToVector(B, a, n , m, pwr + g * a->p + a->t, g * a->p + a->t, h, w, h1);
            }
        }

        pthread_barrier_wait(a->barrier);

        /*if(a->t == 0)
        {
            printf("step %d\n\n\n\n", s);
            
            PrintMatrix(A, a->r, n, n);
            printf("\n\n\n\n");
        }*/
        //pthread_barrier_wait(a->barrier);
        
        for(dg = 2, dob = 1; dg <= pwr; dg *= 2)
        {
            for(i1 = 0; i1 < pwr / dg; i1++)
            {
                //printf("thread %d, step %d, i1 * dg = %d\n", a->t, s, i1 * dg);
                if(a->t - y >= 0 && a->t - y == i1 * dg && x * a->p + a->t + dob < t)
                {
                    //printf("Зануляем (%d, %d)\n", x * a->p + a->t + dob, s);
                    
                    h = (x * a->p + a->t < k) ? m : l;
                    w = (s < k) ? m : l;
                    h1 = (x * a->p + a->t + dob < k) ? m : l;
                
                    if(h1 == w)
                    {
                        PartiallyNullDown(A, a, n, m, x * a->p + a->t + dob, s, x * a->p + a->t, s, h, w, h1);

                        for(z = s + 1; z < t; z++)
                        {
                            w1 = z < k ? m : l;
                            TwicePartiallyApplytoRight(A, a, n, m, x * a->p + a->t + dob, z, x * a->p + a->t, z, h, w, h1, w1);
                        }

                        TwicePartiallyApplyToVector(B, a, n, m, x * a->p + a->t + dob, x * a->p + a->t, h, w, h1);
                    }
                    else
                    {
                        FullNullDown(A, a, n, m, x * a->p + a->t + dob, s, x * a->p + a->t, s, h, w, h1);

                        for(z = s + 1; z < t; z++)
                        {
                            w1 = z < k ? m : l;
                            TwiceApplytoRight(A, a, n, m, x * a->p + a->t + dob, z, x * a->p + a->t, z, h, w, h1, w1);
                        }

                        TwiceApplyToVector(B, a, n, m, x * a->p + a->t + dob, x * a->p + a->t, h, w, h1);
                    }
                }
                if(a->t - y < 0 && a->t - y + a->p == i1 * dg && (x + 1) * a->p + a->t + dob < t)
                {
                    //printf("Зануляем (%d, %d)\n", (x + 1) * a->p + a->t + dob, s);

                    h = (x + 1) * a->p + a->t < k ? m : l;
                    w = (s < k) ? m : l;
                    h1 = (x + 1) * a->p + a->t + dob < k ? m : l;
                
                    if(h1 == w)
                    {
                        PartiallyNullDown(A, a, n, m, (x + 1) * a->p + a->t + dob, s, (x + 1) * a->p + a->t, s, h, w, h1);

                        for(z = s + 1; z < t; z++)
                        {
                            w1 = z < k ? m : l;
                            TwicePartiallyApplytoRight(A, a, n, m, (x + 1) * a->p + a->t + dob, z, (x + 1) * a->p + a->t, z, h, w, h1, w1);
                        }

                        TwicePartiallyApplyToVector(B, a, n, m, (x + 1) * a->p + a->t + dob, (x + 1) * a->p + a->t, h, w, h1);
                    }
                    else
                    {
                        FullNullDown(A, a, n, m, (x + 1) * a->p + a->t + dob, s, (x + 1) * a->p + a->t, s, h, w, h1);

                        for(z = s + 1; z < t; z++)
                        {
                            w1 = z < k ? m : l;
                            TwiceApplytoRight(A, a, n, m, (x + 1) * a->p + a->t + dob, z, (x + 1) * a->p + a->t, z, h, w, h1, w1);
                        }

                        TwiceApplyToVector(B, a, n, m, (x + 1) * a->p + a->t + dob, (x + 1) * a->p + a->t, h, w, h1);
                    }
                }
            }
            dob *= 2;
            pthread_barrier_wait(a->barrier);
        }
    
        //if(a->t == 0)
          //  printf("step = %d\n\n", s);

        //pthread_barrier_wait(a->barrier); 
    }

    //pthread_barrier_wait(a->barrier); 
    /*if(a->t == 0)
    {
        printf("\n\n\n\n");
        PrintMatrix(A, a->r, n, n);
        printf("\n\n\n\n");
    }

    pthread_barrier_wait(a->barrier); */

    for(q = 0; q < hp; q++)
    {
        GetBlock(A, a->block, n, m, q * a->p + a->t, q * a->p + a->t);
        h = (q * a->p + a->t < k) ? m : l;
        for(i = 0, j = 0; i < h; i++)
        {
            c1 = a->block[i * m + i];
            if (std::abs(c1) < a->norm * 1e-15)
            {
                //printf("%e < %e in thread %d, element (%d, %d) of block (%d, %d)\n", c1, a->norm * 1e-15, a->t, i, i, q * a->p + a->t < k, q * a->p + a->t < k);
                a->sol_status = 1;
                j = 1;
                break;
            }
        }

        if(j == 1)
            break;
    }

    if(reduce_sum(a->p, &a->sol_status, 1) > 0)
    {
        a->sol_cpu_time = get_cpu_time() - a->sol_cpu_time;
        a->sol_full_time = get_full_time() - a->sol_full_time;
        a->r1 = -1;
        a->r2 = -1;
        return nullptr;
    }

    //сейчас будет огромный костыль

    if(a->t == 0)
    {
        for(s = t - 1; s >= 0; s--)
        {
            if(t == k + 1)
            {
                if(s == k)
                {
                    FillInv(A, a->T, n, m, s, s);
                    CalcInv(a->T, n, m, s, s);
                    //PrintMatrix(T, 2 * m, m, m * 2);
                    for(i = 0; i < l; i++)
                    {
                        c1 = 0;
                        for(q = 0; q < l; q++)
                            c1 += a->T[i * 2 * m + q + m] * B[s * m + q];
                        X[s * m + i] = c1;
                    }
                }
                else 
                {
                    GetBlock(A, a->block, n, m, s, k);

                    for(i = 0; i < m; i++)
                    {
                        c1 = 0;
                        for(q = 0; q < l; q++)
                            c1 += a->block[i * m + q] * X[k * m + q];
                        a->vec[i] = -c1;
                        a->vec[i] += B[s * m + i];
                    }

                    for(j = k - 1; j > s; j--)
                    {
                        GetBlock(A, a->block, n, m, s, j);

                        for(i = 0; i < m; i++)
                        {
                            c1 = 0;
                            for(q = 0; q < m; q++)
                                c1 += a->block[i * m + q] * X[j * m + q];
                            a->vec[i] -= c1;
                        }
                    }
                    FillInv(A, a->T, n, m, s, s);
                    CalcInv(a->T, n, m, s, s);
                    
                    for(i = 0; i < m; i++)
                    {
                        c1 = 0;
                        for(q = 0; q < m; q++)
                            c1 += a->T[i * 2 * m + q + m] * a->vec[q];
                        X[s * m + i] = c1;
                    }
                }
            }

            if(t == k)
            {
                for(i = 0; i < m; i++)
                    a->vec[i] = B[s * m + i];
                
                for(j = k - 1; j > s; j--)
                {
                    GetBlock(A, a->block, n, m, s, j);

                    for(i = 0; i < m; i++)
                    {
                        c1 = 0;
                        for(q = 0; q < m; q++)
                            c1 += a->block[i * m + q] * X[j * m + q];
                        a->vec[i] -= c1;
                    }
                }

                FillInv(A, a->T, n, m, s, s);
                CalcInv(a->T, n, m, s, s);
                
                for(i = 0; i < m; i++)
                {
                    c1 = 0;
                    for(q = 0; q < m; q++)
                        c1 += a->T[i * 2 * m + q + m] * a->vec[q];
                    X[s * m + i] = c1;
                }
            }
        }
    }
    
    /*if(a->t == 0)
    {
        printf("\n\n\n\n");
        PrintMatrix(A, a->r, n, n);
        printf("\n\n\n\n");
    }*/
    
    a->sol_cpu_time = get_cpu_time() - a->sol_cpu_time;
    pthread_barrier_wait(a->barrier);

    a->sol_full_time = get_full_time() - a->sol_full_time;

    if(a->t == 0)
    {
        printf("Solution vector x:\n");
        PrintMatrix(X, a->r, 1, n);

        if (a->s == 0 && a->filename != nullptr) 
        {
            if(ReadFromFile(a->filename, A, n) != 0)
            {
                a->sol_status = 1;
            }
        }
        else
            InitializeMatrix(A, a->s, n);
        
        CalcB(A, B, n);

        //PrintMatrix(A, n, n, n);
    }
    
    pthread_barrier_wait(a->barrier);

    a->d_cpu_time = get_cpu_time();
    a->d_full_time = get_full_time();


    for(q = 0; q < hp; q++)
    {
        h = (q * a->p + a->t < k) ? m : l;
        for(i = 0; i < h; i++)
            a->r1_norm2 += std::abs(B[(q * a->p + a->t) * m + i]);
    }

    for(q = 0; q < hp; q++)
    {
        for(v = 0; v < t; v++)
        {
            GetBlock(A, a->block, n, m, q * a->p + a->t, v);
            h = (q * a->p + a->t < k) ? m : l;
            w = (v < k) ? m : l;
            for(i = 0; i < h; i++)
            {
                c1 = 0;
                for(j = 0; j < w; j++)
                    c1 += a->block[i * m + j] * X[v * m + j];
                B[(q * a->p + a->t) * m + i] -= c1;
            }
        }
    }

    for(q = 0; q < hp; q++)
    {
        h = (q * a->p + a->t < k) ? m : l;
        for(i = 0; i < h; i++)
            a->r1_norm1 += std::abs(B[(q * a->p + a->t) * m + i]);
    }



    for(q = 0; q < hp; q++)
    {
        h = (q * a->p + a->t < k) ? m : l;
        for(i = 0; i < h; i++)
        {
            c1 = (q * a->p + a->t + i + 1) % 2;
            a->r2_norm1 += std::abs(X[q * a->p + a->t + i] - c1);
            a->r2_norm2 += std::abs(c1);
        }
    }
    
    //printf("before reduce_sum: %.15lf %.15lf\n", a->r2_norm1, a->r2_norm2);

    special_reduce_sum(a->p, &a->r1_norm1, &a->r1_norm2, &a->r2_norm1, &a->r2_norm2, 1);

    special_reduce_sum(a->p, &a->t1, &a->t2, &a->t3, &a->t4, 1);

    a->d_cpu_time = get_cpu_time() - a->d_cpu_time;
    a->d_full_time = get_full_time() - a->d_full_time;
    
    //printf("after reduce_sum: %.20lf %.20lf\n", a->r1_norm1, a->r1_norm2);
    //printf("after reduce_sum: %.20lf %.20lf\n", a->r2_norm1, a->r2_norm2);

    if(a->r1_norm2 > 1e-15 * a->norm)
        a->r1 = a->r1_norm1 / a->r1_norm2;

    if(a->r2_norm2 > 1e-15 * a->norm)
        a->r2 = a->r2_norm1 / a->r2_norm2;
    
    return nullptr;
}