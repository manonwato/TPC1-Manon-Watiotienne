
#include "iterative.h"
#include "jacobi_static.h"
#include "jacobi_dynamic.h"
#include "GaussSeidel_static.h"
#include "GaussSeidel_dynamic.h"
#include "SOR_static.h"
#include "SOR_dynamic.h"
#include "utilities.h"
#include "blas1_static.h"
#include "blas2_static.h"
#include <stdlib.h>


int main (int argc, char** argv){
    
    //=============== déclaration ==============//

    //statiques
    double a_static[Nx][Ny][5]; 
    double b_static[Nx][Nx];
    double u_old_static[Nx][Nx];
    double u_new_static[Nx][Ny];
    double residual_stories[Niter_max];
    // dynamiques
    double*** a_dyn;
    double** b_dyn;
    double** u_old_dyn;
    double** u_new_dyn;
    double* residual_dyn;

    double tol = 1e-10;
    double hx = 1./(Nx-1), hy = 1./(Ny-1);
    double w = 1.941; // pour iterations minimales

     
    //=============== initialisation ==============//

    a_dyn       = allocate_dynamic_array_double_3d(Nx, Ny, 5);
    b_dyn       = allocate_dynamic_array_double_2d(Nx, Ny);
    u_old_dyn   = allocate_dynamic_array_double_2d(Nx, Ny);
    u_new_dyn   = allocate_dynamic_array_double_2d(Nx, Ny);
    residual_dyn = malloc(sizeof(double) * Niter_max);

    init_a(a_static, hx, hy);
    init_a_dyn(a_dyn, hx, hy);
    init_b(b_static, hx, hy);
    init_b_dyn(b_dyn, hx, hy);
    
    set_u(u_old_static, 0.);
    set_u_dyn(u_old_dyn, 0.);
    
    set_u(u_new_static, 0.);
    set_u_dyn(u_new_dyn, 0.);
    dscal(residual_stories, 0., Niter_max);


    // print2D(u_new_static);
    // print2D(b_static);
    // print2D(u_old_static);
    
    // print_a(a_static);
    // print1D(residual_stories, Niter_max);

    //=============== methodes ==============//

    jacobi_static(a_static, b_static, u_old_static, u_new_static, residual_stories, hx, hy, tol);
    set_u(u_old_static, 0.);  
    set_u(u_new_static, 0.);

    jacobi_dynamic( a_dyn, b_dyn, u_old_dyn, u_new_dyn, residual_dyn, hx, hy, tol);
    set_u_dyn(u_old_dyn, 0.);
    set_u_dyn(u_new_dyn, 0.);

    GaussSeidel_static(a_static, b_static, u_old_static, u_new_static, residual_stories, hx, hy, tol);
    set_u(u_old_static, 0.);  
    set_u(u_new_static, 0.);

    GaussSeidel_dynamic(a_dyn, b_dyn, u_old_dyn, u_new_dyn, residual_dyn, hx, hy, tol);
    set_u_dyn(u_old_dyn, 0.);
    set_u_dyn(u_new_dyn, 0.);

    SOR_static(a_static, b_static, u_old_static, u_new_static, residual_stories, hx, hy, w, tol);
    
    SOR_dynamic(a_dyn, b_dyn, u_old_dyn, u_new_dyn, residual_dyn, hx, hy, w, tol);

    // print1D(residual_stories, Niter_max);
    // print2D(u_new_static);

    //=============== free dynamique ==============//
    free_dynamic_array_double_3d(a_dyn, Nx, Ny, 5);
    free_dynamic_array_double_2d(b_dyn, Nx, Ny);
    free_dynamic_array_double_2d(u_old_dyn, Nx, Ny);
    free_dynamic_array_double_2d(u_new_dyn, Nx, Ny);
    free(residual_dyn);
    return 0;
}