#include "header.h"

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

double CalcR1(double* A, const double* X, double* B, double* block, int n, int m) 
{
    double norm1 = 0, norm2 = 0;
    double tmp;

    int k = n / m;
    int l = n % m;

    int i, j, p, q, h, w;
    int t = (l == 0) ? k : k + 1; 

    for(i = 0; i < n; i++)
        norm2 += std::abs(B[i]);

    for(i = 0; i < t; i++)
    {
        for(j = 0; j < t; j++)
        {
            GetBlock(A, block, n, m, i, j);
            h = (i < k) ? m : l;
            w = (j < k) ? m : l;
            for(p = 0; p < h; p++)
            {
                tmp = 0;
                for(q = 0; q < w; q++)
                    tmp += block[p * m + q] * X[j * m + q];
                B[i * m + p] -= tmp;
            }
        }
    }

    for(i = 0; i < n; i++)
        norm1 += std::abs(B[i]);
    
    return norm1 / norm2;
}

double CalcR2(const double* X, int n) 
{
    double norm1 = 0, norm2 = 0;
    double x1;
    for (int i = 0; i < n; i++) 
    {
        x1 = double((i + 1) % 2);
        norm1 += std::abs(X[i] - x1);
        norm2 += std::abs(x1);
    }
    return norm1 / norm2;
}

int SaulGoodman(double* A, double *B, double* X, double* block, double* T, double* mat, double* vec, int n, int m)
{
    int k = n / m;
    int l = n % m;

    int t = (l == 0) ? k : k + 1;

    //printf("t = %d\n", t);

    int h, w, s, i, j, q, v, z, cnt, cnt1;

    double cos, sin, sq, c1, c2, eps; 

    double norm = A[0];
    
    for(int i = 1; i < n; i++)
    {
        c1 = A[i];
        if(norm > c1)
            norm = c1;
    }

    eps = 1e-15 * norm;

    //PrintMatrix(B, n, 1, n);

    //double tmp;
    //tmp = clock();

    for(s = 0; s < t; s++)
    {
        //tmp = clock();
        h = w = (s < k) ? m : l;
        //printf("s = %d\n", s);
        GetBlock(A, block, n, m, s, s);
        //PrintMatrix(block, h, h, h);
        for(j = 0; j < w; j++)
        {
            for(i = j + 1; i < h; i++)
            {
                cnt1++;
                c1 = block[j * m + j];
                c2 = block[i * m + j];
                
                sq = sqrt(c1 * c1  + c2 * c2);
                //printf("i = %d j = %d c1 = %10.3e c2 = %10.3e\n", i, j, c1, c2);
                //PrintMatrix(block, h, h, h);
                if(sq > eps)
                {
                    cnt++;
                    cos = c1 / sq;
                    sin = -c2 / sq;
                    T[i * 2 * m + j] = cos;
                    T[i * 2 * m + j + m] = sin;

                    block[i * m + j] = 0;
                    block[j * m + j] = sq;
                    for(q = j + 1; q < w; q++) // j + 1
                    {
                        c1 = block[i * m + q]; 
                        c2 = block[j * m + q]; 
                        block[i * m + q] = cos * c1 + sin * c2;
                        block[j * m + q] = cos * c2 - sin * c1;
                    }
                }
                else
                    return -1;
            }
            
        }
        //if(std::abs(block[(w - 1) * m + w -1]) < eps)
          //  return -1;
            
        SetBlock(A, block, n, m, s, s);

        //printf("s = %d: %lf diag->diag\n", s, (clock() - tmp)/CLOCKS_PER_SEC);
        //tmp = clock();
        for(v = s + 1; v < t; v++)
        {
            GetBlock(A, block, n, m, s, v);
            w = (v < k) ? m : l;
            for(j = 0; j < m; j++)
            {
                for(i = j + 1; i < m; i++)
                {
                    for(q = 0; q < w; q++)
                    {
                        c1 = block[i * m + q];
                        c2 = block[j * m + q];
                        cos = T[i * 2 * m + j]; 
                        sin = T[i * 2 * m + j + m];
                        block[i * m + q] = c1 * cos + c2 * sin;
                        block[j * m + q] = c2 * cos - c1 * sin;
                    }
                }
            }
            SetBlock(A, block, n, m, s, v);
        }

        //printf("s = %d: %lf to right blocks\n", s,  (clock() - tmp)/CLOCKS_PER_SEC);
        //tmp = clock();
        h = (s < k) ? m : l;
        //printf("s = %d\n", s);
        for(j = 0; j < h; j++)
        {
            for(i = j + 1; i < h; i++)
            {
                c1 = B[i + s * m];
                c2 = B[j + s * m];
                //printf("i = %d j = %d c1 = %10.3e c2 = %10.3e\n", i, j, c1, c2);
                cos = T[i * 2 * m + j];
                sin = T[i * 2 * m + j + m];
                B[i + s * m] = cos * c1 + sin * c2;
                B[j + s * m] = cos * c2 - sin * c1;
            }
        }
        //PrintMatrix(A, n, n, n);
        //printf("\n");

        //printf("s = %d: %lf to rigth blocks and row\n", s, (clock() - tmp)/CLOCKS_PER_SEC);
        //tmp = clock();
        for(v = s + 1; v < t; v++)
        {
            //tmp = clock();
            h = (v < k) ? m : l;
            aboba(A, mat, n, m, v, s, s, s);
            for(j = 0; j < m; j++)  
            {
                for(i = 0; i < h; i++)
                {
                    c1 = mat[j * m + j];
                    c2 = mat[i * m + j + m * m];
                    //printf("v = %d s = %d i = %d j = %d c1 = %10.3e c2 = %10.3e\n", v, s, i, j, c1, c2);
                    //PrintMatrix(mat, h + 1, h + 1, m);
                    sq = sqrt(c1 * c1 + c2 * c2);
                    if(sq > eps)
                    {
                        cos = c1 / sq;
                        sin = -c2 / sq;
                        T[i * 2 * m + j] = cos;
                        T[i * 2 * m + j + m] = sin;

                        mat[i * m + j + m * m] = 0;
                        mat[j * m + j] = sq;

                        for(q = j + 1; q < m; q++) 
                        {
                            c1 = mat[j * m + q];
                            c2 = mat[i * m + q + m * m];    
                            mat[j * m + q] = cos * c1 - sin * c2;
                            mat[i * m + q + m * m] = cos * c2 + sin * c1;
                        }
                    }
                }
            }
            aboba_set(A, mat, n, m, v, s, s, s);
            //PrintMatrix(A, n, n, n);
            //printf("\n");
            //printf("%lf\n", (clock() - tmp)/CLOCKS_PER_SEC);
            //tmp = clock();
            
            for(z = s + 1; z < t; z++)
            {
                w = z < k ? m : l;
                aboba(A, mat, n, m, v, z, s, z);
                for(j = 0; j < m; j++)
                {
                    for(i = 0; i < h; i++)
                    {
                        cos = T[i * 2 * m + j];
                        sin = T[i * 2 * m + j + m];
                        for(q = 0; q < w; q++)
                        {
                            c1 = mat[j * m + q];
                            c2 = mat[i * m + q + m * m];
                            mat[j * m + q] = cos * c1 - sin * c2;
                            mat[i * m + q + m * m] = cos * c2 + sin * c1; 
                        }
                    }
                }
                aboba_set(A, mat, n, m, v, z, s, z);
            }

            //PrintMatrix(A, n, n, n);
            //printf("\n");

            //printf("%lf\n", (clock() - tmp)/CLOCKS_PER_SEC);
            //tmp = clock();
            aboba1(B, vec, n, m, v, s);
            for(j = 0; j < m; j++)
            {
                for(i = 0; i < h; i++)
                {
                    c1 = vec[j];
                    c2 = vec[i + m];

                    cos = T[i * 2 * m + j];
                    sin = T[i * 2 * m + j + m];

                    vec[j] = cos * c1 - sin * c2;
                    vec[i + m] = cos * c2 + sin * c1;
                }
            }
            aboba1_set(B, vec, n, m, v, s);
            //printf("%lf\n", (clock() - tmp)/CLOCKS_PER_SEC);
        }
        //printf("s = %d: %lf - to null below\n", s, (clock() - tmp)/CLOCKS_PER_SEC);
    }
    //PrintMatrix(A, n, n, n);
    //std::cout << std::endl;
    //PrintMatrix(B, n, 1, n);

    //printf("%lf\n", (clock() - tmp)/CLOCKS_PER_SEC);

    //tmp = clock();

    for(s = t - 1; s >= 0; s--)
    {
        if(t == k + 1)
        {
            if(s == k)
            {
                FillInv(A, T, n, m, s, s);
                CalcInv(T, n, m, s, s);
                //PrintMatrix(T, 2 * m, m, m * 2);
                for(i = 0; i < l; i++)
                {
                    c1 = 0;
                    for(q = 0; q < l; q++)
                        c1 += T[i * 2 * m + q + m] * B[s * m + q];
                    X[s * m + i] = c1;
                }
            }
            else 
            {
                GetBlock(A, block, n, m, s, k);

                for(i = 0; i < m; i++)
                {
                    c1 = 0;
                    for(q = 0; q < l; q++)
                        c1 += block[i * m + q] * X[k * m + q];
                    vec[i] = -c1;
                    vec[i] += B[s * m + i];
                }

                for(j = k - 1; j > s; j--)
                {
                    GetBlock(A, block, n, m, s, j);

                    for(i = 0; i < m; i++)
                    {
                        c1 = 0;
                        for(q = 0; q < m; q++)
                            c1 += block[i * m + q] * X[j * m + q];
                        vec[i] -= c1;
                    }
                }
                FillInv(A, T, n, m, s, s);
                CalcInv(T, n, m, s, s);
                
                for(i = 0; i < m; i++)
                {
                    c1 = 0;
                    for(q = 0; q < m; q++)
                        c1 += T[i * 2 * m + q + m] * vec[q];
                    X[s * m + i] = c1;
                }
            }
        }

        if(t == k)
        {
            for(i = 0; i < m; i++)
                vec[i] = B[s * m + i];
            
            for(j = k - 1; j > s; j--)
            {
                GetBlock(A, block, n, m, s, j);

                for(i = 0; i < m; i++)
                {
                    c1 = 0;
                    for(q = 0; q < m; q++)
                        c1 += block[i * m + q] * X[j * m + q];
                    vec[i] -= c1;
                }
            }

            FillInv(A, T, n, m, s, s);
            CalcInv(T, n, m, s, s);
            
            for(i = 0; i < m; i++)
            {
                c1 = 0;
                for(q = 0; q < m; q++)
                    c1 += T[i * 2 * m + q + m] * vec[q];
                X[s * m + i] = c1;
            }
        }
    }

    //printf("%lf\n", (clock() - tmp)/CLOCKS_PER_SEC);
    //printf("\n\n\n\n");
    //PrintMatrix(A, n, n, n); 
    //printf("\n\n\n\n");
    return 0;
}