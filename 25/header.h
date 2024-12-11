#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <cmath>
#include <sys/time.h>
#include <sys/resource.h>

#ifndef HEADER_H
#define HEADER_H

int ReadFromFile(const char* filename, double* matrix, int n);
void InitializeMatrix(double* A, int k, int n);
void PrintMatrix(const double* A, int m, int l, int n); 
double get_full_time(void);
void GetInv(double* A, double* norm, double* trace, double* len, int n);
bool IsSymmetric(double* A, double norm, double eps, int n);
int GetEV(double* A, double* X, double* A1, double* Q, bool* mask, double norm, double eps, int n, int* it);

void Make3Diag(double* A, double norm, double eps, int n);
void BuildQR(double* A, double* Q, bool* mask, bool* flag1, double norm, double eps, int n, int k);
void FindNext(double* A, double* Q, bool* mask, double norm, double eps, int n, int k);

double CalcR1(double* X, double trace, double norm, int n);
double CalcR2(double* X, double len, double norm, int n);

#endif 