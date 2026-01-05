#pragma once
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "iterative.h"

void GaussSeidel_dynamic(double*** a, double** b, double** u_old, double** u_new, double* residual, double hx, double hy, double tol);
