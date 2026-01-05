#pragma once
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "iterative.h"

void SOR_static(
        double a[Nx][Ny][5], 
        double b[Nx][Ny], 
        double u_old[Nx][Ny], 
        double u_new[Nx][Ny], 
        double residual[Niter_max], 
        double hx, 
        double hy, 
        double w, 
        double tol);
