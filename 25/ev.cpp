#include "header.h"

#define MAX_ITERS 20000

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

void InitializeMatrix(double* A, int k, int n)
{
    auto f = [](int k, int i, int j, int n) -> double
    {
        switch(k)
        {
            case 1: return n - std::max(i, j);
            case 2: if(i == j) return 2; if(std::abs(i - j) == 1) return -1; return 0;
            case 3: if(i == j && j < n - 1) return 1; if(j == n - 1) return i + 1; if(i == n - 1) return j + 1; return 0;
            case 4: return 1.0 / (i + j + 1);
            //case 5: if(i == j) return i; else return 0;
            default: return 0;
        }
    };

    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            A[i * n + j] = f(k, i, j, n);
}

void PrintMatrix(const double* A, int m, int l, int n) 
{
    for (int i = 0; i < std::min(m, l); ++i) 
    {
        for (int j = 0; j < std::min(m, n); ++j) 
        {
            printf(" %10.3e", A[i * n + j]);
        }

        printf("\n");
    }

    printf("\n");
}

double get_full_time() 
{
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec/1e6;
}

double CalcR1(double* X, double trace, double norm, int n)
{
    int i;
    double s = 0;

    for(i = 0; i < n; i++)
        s += X[i];

    return std::abs(trace - s) / norm;
}

double CalcR2(double* X, double len, double norm, int n)
{
    int i;
    double s = 0;

    for(i = 0; i < n; i++)
        s += X[i] * X[i];


    return std::abs(len - sqrt(s)) / norm;
}

void GetInv(double* A, double* norm, double* trace, double* len, int n)
{
    int i, j;
    double tmp = 0, tr = 0, len1 = 0, nrm = 0;
    
    for(i = 0; i < n; i++)
        nrm += std::abs(A[i]);
    
    for(i = 1; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            tmp += std::abs(A[i * n + j]);
        }
        if(tmp > nrm)
            nrm = tmp;
        tmp = 0;
    }

    for(i = 0; i < n; i++)
        tr += A[i * n + i];

    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            len1 += A[i * n + j] * A[j * n + i]; 

    len1 = sqrt(len1);

    *norm = nrm;
    *trace = tr;
    *len = len1;
}

bool IsSymmetric(double* A, double norm, double eps, int n)
{
    int i, j;

    double normeps = norm * eps;

    for(i = 0; i < n; i++)
    {
        for(j = i + 1; j < n; j++)
        {
            if(std::abs(A[i * n + j] - A[j * n + i]) > normeps)
                return false;
        }
    }
    return true;
}

void Make3Diag(double* A, double norm, double eps, int n)
{
    int i, j, q;

    double normeps = norm * eps;

    double cos, sin, sq, c1, c2, c3, c4;

    for(j = 0; j < n; j++)
    {
        for(i = j + 2; i < n; i++)
        {
            c1 = A[(j + 1) * n + j];
            c2 = A[i * n + j];
            
            sq = sqrt(c1 * c1  + c2 * c2);
            //printf("i = %d j = %d c1 = %10.3e c2 = %10.3e\n", i, j, c1, c2);
            //PrintMatrix(block, h, h, h);
            if(sq > normeps)
            {
                cos = c1 / sq;
                sin = -c2 / sq;

                A[i * n + j] = 0;
                A[(j + 1) * n + j] = sq;
                for(q = j + 1; q < n; q++) // j + 1
                {
                    c1 = A[i * n + q]; 
                    c2 = A[(j + 1) * n + q]; 
                    A[i * n + q] = cos * c1 + sin * c2;
                    A[(j + 1) * n + q] = cos * c2 - sin * c1;
                }
                
                //printf("%d %d\n", i, j + 1);

                //printf("before mult from right:\n");
                //PrintMatrix(A, n, n, n);
                c1 = A[i * n + i];
                c2 = A[i * n + j + 1];
                c3 = A[(j + 1) * n + i];
                c4 = A[(j + 1) * n + j + 1];


                A[i * n + i] = c1 * cos + c2 * sin;
                A[i * n + j + 1] = -(c1 * sin - c2 * cos);
                A[(j + 1) * n + i] = c3 * cos + c4 * sin;
                A[(j + 1) * n + j + 1] = -(c3 * sin - c4 * cos);
                

                A[j * n + i] = 0;
                A[j * n + j + 1] = sq;

                for(q = j + 1; q < n; q++)
                {
                    if(q != i && q != j + 1)
                    {
                        A[q * n + i] = A[i * n + q];
                        A[q * n + j + 1] = A[(j + 1) * n + q];
                    }
                }
                //printf("(%d, %d) = 0 (%d, %d) = sq\n", j, i, j, j + 1);
                /*for(q = j + 1; q < n; q++)
                {
                    c1 = A[q * n + i]; 
                    c2 = A[q * n + j + 1]; 
                    A[q * n + i] = cos * c1 + sin * c2;
                    A[q * n + j + 1] = cos * c2 - sin * c1;
                    
                }*/

                //printf("after mult from right:\n");
                //PrintMatrix(A, n, n, n);
            }
            else
            {
                A[i * n + j] = 0;
                A[(j + 1) * n + j] = sq;
                A[j * n + i] = 0;
                A[j * n + j + 1] = sq;
            }

            
        }
    }
}

void BuildQR(double* A, double* Q, bool* mask, bool* flag1, double norm, double eps, int n, int k)
{
    int i, j, q;

    double normeps = norm * eps;

    double cos, sin, sq, c1, c2;

    //printf("im in QR at k = %d\n", k);

    for(j = 0; j < k - 1; j++)
    {
        for(i = j + 1; i < j + 2; i++)
        {
            c1 = A[j * n + j];
            c2 = A[i * n + j];
            
            sq = sqrt(c1 * c1  + c2 * c2);
            //printf("i = %d j = %d c1 = %10.3e c2 = %10.3e\n", i, j, c1, c2);
            //PrintMatrix(block, h, h, h);
            if(sq > normeps)
            {
                mask[j] = true;
                cos = c1 / sq;
                sin = -c2 / sq;
                Q[j] = cos;
                Q[n + j] = -sin;

                A[i * n + j] = 0;
                A[j * n + j] = sq;
                for(q = j + 1; q < std::min(j + 3, k); q++) // j + 1
                {
                    c1 = A[i * n + q]; 
                    c2 = A[j * n + q]; 
                    A[i * n + q] = cos * c1 + sin * c2;
                    A[j * n + q] = cos * c2 - sin * c1;
                }
            }
            else
            {
                mask[j] = false;
                A[i * n + j] = 0;
                A[j * n + j] = sq;
                Q[j] = 0;
                Q[n + j] = 0;
                
                *flag1 = true;
                return;
            }
        }
    }
}

void FindNext(double* A, double* Q, bool* mask, double norm, double eps, int n, int k)
{
    int i, j, q;

    double normeps = norm * eps;

    double cos, sin, c1, c2, c3, c4;

    for(j = 0; j < k - 1; j++)
    {
        for(i = j + 1; i < j + 2; i++)
        {
            if(mask[j] == true)
            {
                for(q = std::max(0, j - 1); q < std::min(k, j + 2); q++)
                {
                    c1 = A[q * n + i];
                    c2 = A[q * n + j];
                    cos = Q[j]; 
                    sin = -Q[n + j];
                    c3 = c1 * cos + c2 * sin;
                    c4 = c2 * cos - c1 * sin;
                    if(std::abs(c3) > normeps) A[q * n + i] = c3; else A[q * n + i] = 0;
                    if(std::abs(c4) > normeps) A[q * n + j] = c4; else A[q * n + j] = 0;
                }
            }
        }
    }
}


int GetEV(double* A, double* X, double* A1, double* Q, bool* mask, double norm, double eps, int n, int* it)
    {
    int iters = 0, i;

    int k = n;

    double s, normeps = norm * eps;

    //Make3Diag(A, norm, eps, n);

    //bool flag = false;

    bool flag1 = false;


    while(iters < MAX_ITERS)
    {
        if(k == 1)
        {
            X[0] = A[0];
            break;
        }
        if(k == 2)
        {
            s = sqrt((A[0] - A[n + 1]) * (A[0] - A[n + 1]) + 4 * A[1] * A[1]);
            X[0] = (A[0] + A[n + 1] + s) / 2;
            X[1] = (A[0] + A[n + 1] - s) / 2;
            break;
        }
        if(std::abs(A[(k - 1) * n + k - 2]) < normeps)
        {
            //printf("GOT EV HERE\n nulling (%d, %d)\n", k - 1, k - 2);
            //PrintMatrix(A, k, k, k);
            A[(k - 1) * n + k - 2] = 0;
            X[k - 1] = A[(k - 1) * n + k - 1];
            k--;
            iters++;
            continue;
        }

        for(i = 0; i < k - 1; i++)
            A1[i] = A[i * n + i + 1];
    
        for(i = 0; i < k; i++)
            A1[n + i] = A[i * n + i];

        for(i = 1; i < k; i++)
            A1[2 * n + i - 1] = A[i * n + i - 1];
        
        //if(!flag)
            s = A[(k - 1) * n + k - 1] + A[(k - 1) * n + k - 2] / 2;
        /*else
        {
            s = A[(k - 1) * n + k - 1] - A[(k - 1) * n + k - 2] / 2;
            flag = false;
        }*/

        //s = 0.283119;
        //s = A[(k - 1) * n + k - 1];

        /*trace = 0;
        for(i = 0; i < n; i++)
            trace += A[i * n + i];*/

        //printf("k: %d, s = %lf, trace = %lf, it = %d\n", k, s, trace, iters);

        //PrintMatrix(A, n, n, n);

        for(i = 0; i < k; i++)
            A[i * n + i] -= s;

        //printf("Before QR\n");
        //PrintMatrix(A, n, n, n);
        BuildQR(A, Q, mask, &flag1, norm, eps, n, k);

        //printf("After QR\n");
        //PrintMatrix(A, n, n, n);
        
        if(flag1)
        {
            //printf("ugadali na k = %d, s = %lf, it = %d\n", k, s, iters);

            for(i = 0; i < k - 1; i++)
                A[i * n + i + 1] = A1[i];
    
            for(i = 0; i < k; i++)
                A[i * n + i] = A1[n + i];

            for(i = 1; i < k; i++)
                A[i * n + i - 1] = A1[2 * n + i - 1];

            X[k - 1] = s;
            k--;
            iters++;
            //flag = true;
            flag1 = false;
            continue;
        }
        FindNext(A, Q, mask, norm, eps, n, k);

        for(i = 0; i < k; i++)
            A[i * n + i] += s;


        for(i = 0; i < k - 1; i++)
            A1[i] = A[i * n + i + 1];
    
        for(i = 0; i < k; i++)
            A1[n + i] = A[i * n + i];

        for(i = 1; i < k; i++)
            A1[2 * n + i - 1] = A[i * n + i - 1];

        iters++;
    }

    //printf("number of iterations: %d\n", iters);

    *it = iters;

    if(iters == MAX_ITERS)
        return -1;
    return 0;
}

