#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <cmath>

#ifndef HEADER_H
#define HEADER_H

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
double CalcR1(double* A, const double* X, double* B, double* block, int n, int m);
double CalcR2(const double* X, int n);
int SaulGoodman(double* A, double *B, double* X, double* block, double* T, double* mat, double* vec, int n, int m);

#endif 