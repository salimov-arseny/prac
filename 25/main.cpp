#include "header.h"

#include <fenv.h>

int main(int argc, char* argv[]) 
{
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    if (argc < 5) {
        printf("Usage: %s n m eps k [filename]\n", argv[0]);
        return 1;
    }

    int n, m, k, iters = 0;
    double eps;
    double* A = nullptr, *X = nullptr, *A1 = nullptr, *Q = nullptr;
    bool* mask = nullptr;

    try
    {
        n = std::stoi(argv[1]);
        m = std::stoi(argv[2]);
        eps = std::stod(argv[3]);
        k = std::stoi(argv[4]);
    }
    catch(const std::invalid_argument&)
    {
        printf("Wrong arguments\n"); 
        return 1;
    }

    if(n <= 0 || m < 0 || k < 0 || k > 4)
    {
        printf("Wrong arguments\n"); 
        return 1;
    }

    int sol_status;

    double r1, r2, t1, t2, norm, trace, len;

    try
    {
        A = new double[n * n];
        X = new double[n];
        A1 = new double[3 * n];
        Q = new double[2 * n];
        mask = new bool[n];
    }
    catch(const std::bad_alloc&)
    {
        printf("Not enough memory\n");
        delete[] A;
        delete[] X;
        delete[] A1;
        delete[] Q;
        delete[] mask;
        return 1;
    }
    

    memset(A1, 0, 3 * n * sizeof(double));
    memset(Q, 0, 2 * n * sizeof(double));
    memset(mask, 0, n * sizeof(bool));


    if (k == 0 && argc == 6) 
    {
        const char* filename = argv[5];
        if(ReadFromFile(filename, A, n) != 0)
        {
            delete[] A;
            delete[] X;
            delete[] A1;
            delete[] Q;
            delete[] mask;
            return 2;
        }
    }
    else
        InitializeMatrix(A, k, n);

    printf("Matrix A:\n");

    PrintMatrix(A, m, n, n);

    GetInv(A, &norm, &trace, &len, n);

    //printf("norm = %lf trace = %lf len = %lf\n", norm, trace, len);

    if(!IsSymmetric(A, norm, eps, n))
    {
        printf("A is not symmetric\n");
        printf ("%s : Residual1 = %d Residual2 = %d Iterations = %d Iterations1 = %d Elapsed1 = %d Elapsed2 = %d\n",
        argv[0], -1, -1, iters, iters / n, 0, 0);
        delete[] A;
        delete[] X;
        delete[] A1;
        delete[] Q;
        delete[] mask;
        return 3;
    }

    t1 = get_full_time();
    Make3Diag(A, norm, eps, n);

    //PrintMatrix(A, n, n, n);
    t1 = get_full_time() - t1;

    t2 = get_full_time();

    sol_status = GetEV(A, X, A1, Q, mask, norm, eps, n, &iters);

    t2 = get_full_time() - t2;

    
    if(sol_status != 0) 
    {
        printf("Exceeded MAX_ITERS\n");
        delete[] A;
        delete[] X;
        delete[] A1;
        delete[] Q;
        delete[] mask;

        printf ("%s : Residual1 = %d Residual2 = %d Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
        argv[0], -1, -1, iters, iters / n, t1, t2);

        return 1;
    }

    r1 = CalcR1(X, trace, norm, n);
    r2 = CalcR2(X, len, norm, n);

    printf("Eigenvalues:\n");
    PrintMatrix(X, m, 1, n);

    if (k == 0 && argc == 6) 
    {
        const char* filename = argv[5];
        if(ReadFromFile(filename, A, n) != 0)
        {
            delete[] A;
            delete[] X;
            delete[] A1;
            delete[] Q;
            delete[] mask;
            return 1;
        }
    }
    else
        InitializeMatrix(A, k, n);
    

    printf ("%s : Residual1 = %e Residual2 = %e Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
    argv[0], r1, r2, iters, iters / n, t1, t2);

    delete[] A;
    delete[] X;
    delete[] A1;
    delete[] Q;
    delete[] mask;

    return 0;
}

