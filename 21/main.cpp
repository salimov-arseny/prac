#include "header.h"

int main(int argc, char* argv[]) 
{
    if (argc < 5) {
        printf("Usage: %s n m r s [filename]\n", argv[0]);
        return 1;
    }

    int n = std::atoi(argv[1]); 
    int m = std::atoi(argv[2]); 
    int r = std::atoi(argv[3]); 
    int s = std::atoi(argv[4]); 

    int sol_status;

    double r1, r2;

    std::chrono::duration<double> t1, t2;

    double* A = new double[n * n];
    double* B = new double[n];
    double* X = new double[n]; 
    double* block = new double[m * m];
    double* T = new double[m * 2 * m];
    double* mat = new double[2 * m * m];
    double* vec = new double[2 * m];

    if (s == 0 && argc == 6) 
    {
        const char* filename = argv[5];
        if(ReadFromFile(filename, A, n) != 0)
        {
            delete[] A;
            delete[] B;
            delete[] X;
            delete[] block;
            delete[] T;
            delete[] mat;
            delete[] vec;
            return 1;
        }
    }
    else
        InitializeMatrix(A, s, n);

    printf("Matrix A:\n");

    PrintMatrix(A, r, n, n);

    CalcB(A, B, n);

    auto start = std::chrono::high_resolution_clock::now();

    sol_status = SaulGoodman(A, B, X, block, T, mat, vec, n, m);

    t1 = std::chrono::high_resolution_clock::now() - start;

    if (sol_status != 0) 
    {
        printf("Error: Unable to solve the system\n");
        delete[] A;
        delete[] B;
        delete[] X;
        delete[] block;
        delete[] T;
        delete[] mat;
        delete[] vec;

        printf("%s : Task = %d Res1 = %d Res2 = %d T1 = %.2f T2 = %d S = %d N = %d M = %d\n", argv[0], 21, -1, -1, t1.count(), 0, s, n, m);

        return 1;
    }

    printf("Solution vector x:\n");
    PrintMatrix(X, r, 1, n);

    if (s == 0 && argc == 6) 
    {
        const char* filename = argv[5];
        if(ReadFromFile(filename, A, n) != 0)
        {
            delete[] A;
            delete[] B;
            delete[] X;
            delete[] block;
            delete[] T;
            delete[] mat;
            delete[] vec;
            return 1;
        }
    }
    else
        InitializeMatrix(A, s, n);
    
    CalcB(A, B, n);

    start = std::chrono::high_resolution_clock::now();

    r1 = CalcR1(A, X, B, block, n, m);
    r2 = CalcR2(X, n);

    t2 = std::chrono::high_resolution_clock::now() - start;

    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], 21, r1, r2, t1.count(), t2.count(), s, n, m);

    delete[] A;
    delete[] B;
    delete[] X;
    delete[] block;
    delete[] T;
    delete[] mat;
    delete[] vec;

    return 0;
}

