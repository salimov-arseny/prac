#include "header.h"

#include <fenv.h>

int main(int argc, char* argv[]) 
{
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    if (argc < 6) {
        printf("Usage: %s n m p r s [filename]\n", argv[0]);
        return 1;
    }
    int n, m, p, r, s;

    char* filename;

    try
    {
        n = std::stoi(argv[1]); 
        m = std::stoi(argv[2]);
        p = std::stoi(argv[3]);
        r = std::stoi(argv[4]); 
        s = std::stoi(argv[5]);
        filename = argv[6]; 
    }
    catch(const std::invalid_argument&)
    {
        printf("Wrong arguments\n"); 
        return 1;
    }

    int sol_status;

    double r1, r2;

    double t1, t2;

    double* A = new double[n * n];
    double* B = new double[n];
    double* X = new double[n]; 

    Args* args = new Args[p];

    pthread_barrier_t barrier;

    if (s == 0 && argc == 7) 
        filename = argv[6];
    else
        filename = nullptr;

    pthread_barrier_init(&barrier, 0, p);

    for(int i = 0; i < p; i++)
    {
        if(args[i].init(n, m, r, s, filename, i, p, argc, A, B, X, &barrier) != 0)
        {
            //printf("Not enough memory for thread %d data\n", i);
            delete[] A;
            delete[] B;
            delete[] X;
            delete[] args;
            return 2;
        }
    }

    for(int i = 1; i < p; i++)
    {
        if(pthread_create(&args[i].tid, nullptr, thread_func, args + i))
        {
            printf("Can not create thread %d\n", i);
            delete[] A;
            delete[] B;
            delete[] X;
            delete[] args;
            return 3;
        }
    }

    args[0].tid = pthread_self();
    thread_func(args + 0);
    
    for(int i = 1; i < p; i++)
        pthread_join(args[i].tid, nullptr);

    //PrintMatrix(A, n, n, n);


    sol_status = args[0].sol_status;
    t1 = args[0].sol_full_time;
    t2 = args[0].d_full_time;
    r1 = args[0].r1;
    r2 = args[0].r2;

    if (sol_status > 0) 
    {
        printf("Error: Unable to solve the system\n");
        delete[] A;
        delete[] B;
        delete[] X;
        delete[] args;

        printf("%s : Task = %d Res1 = %d Res2 = %d T1 = %.2f T2 = %d S = %d N = %d M = %d P = %d\n", argv[0], 21, -1, -1, t1, 0, s, n, m, p);

        return 1;
    }
    else if(sol_status < 0)
    {
        printf("Error: Cant open file\n");
        delete[] A;
        delete[] B;
        delete[] X;
        delete[] args;

        printf("%s : Task = %d Res1 = %d Res2 = %d T1 = %.2f T2 = %d S = %d N = %d M = %d P = %d\n", argv[0], 21, -1, -1, t1, 0, s, n, m, p);

        return 1;
    }

    //for(int i = 0; i < p; i++)
      //  args[i].print();

    double s1 = 0;
    double s2 = 0;

    for(int i = 0; i < p; i++)
    {
        s1 += args[i].sol_cpu_time * args[i].sol_cpu_time;
        s2 += args[i].sol_cpu_time;
    }

    //printf("cpu_time variance = %lf\n", s1/p - s2*s2/(p*p));

    s1 = s2 = 0;

    for(int i = 0; i < p; i++)
    {
        s1 += args[i].sol_full_time * args[i].sol_full_time;
        s2 += args[i].sol_full_time;
    }
    //printf("full_time variance = %lf\n", s1/p - s2*s2/(p*p));

    //printf("t1 = %lf\nt2 = %lf\nt3 = %lf\n", args[0].t1, args[0].t2, args[0].t3);

    printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], 21, r1, r2, t1, t2, s, n, m, p);
    
    delete[] A;
    delete[] B;
    delete[] X;
    delete[] args;
    pthread_barrier_destroy(&barrier);

    return 0;
}
