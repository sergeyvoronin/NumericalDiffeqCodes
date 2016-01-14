/* simple multigrid solver 

compile with: g++ -o mgs mg_solver.cpp -O2
run with ./mgs
make sure there is a directory called data in same location as executable

Sergey Voronin 
*/

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

//function prototypes
void startup();
void gauss_seidel_iter(double **U, double **rhs, int level);
void coarseToFine(double **coarse_vec, double **fine_vec, int level);
void fineToCoarse(double **fine_vec, double **coarse_vec, int level);
void mg_vcycle(double **vec, double **rhs, int level, int num_pre_iters, int num_post_iters);
double bndry_func(double x, double y);
double rhs_func(double x, double y);
double sol_func(double x, double y);

//global vars
int *grid_sizes;
double *step_sizes;
int xi, xf, yi, yf, num_grid_levels, problem_num;


int main(){
    // set up vars and arrays ------>
    int i,j, max_size, num_pre_iters, num_post_iters, num_v_cycles, n;
    double x,y,h,true_norm,solution_norm,error_norm,percent_error;
    double **U_mat;
    double **rhs_mat;

    // initialize domain
    xi = 0; xf = 1;
    yi = 0; yf = 1;

    // pick a test problem to use
    // omega = {(x,y): 0<x<1; 0<y<1}
    // problem 1: u_xx + u_yy = -2*pi^2*sin(pi*x)*sin(pi*y) on omega; u(x,y) = 0 on bdry
    // problem 2: u_xx + u_yy = 4 on omega; u(x,y) = x^2 + y^2 pon bdry
    //problem_num = 2; 
    problem_num = 2; 

    // set multigrid params
    num_grid_levels = 5;
    num_v_cycles = 5;
    num_pre_iters = 50;
    num_post_iters = 50;

    printf("performing startup stuff..\n");
    startup();

    // set up grid sizes and arrays
    // make space for arrays ----->
    // solution on fine grid
    printf("allocating space for solution mat..\n");
    max_size = grid_sizes[num_grid_levels-1];
    U_mat = (double**)malloc(max_size*sizeof(double*));
    for(i=0; i<max_size; i++){ 
        U_mat[i] = (double*)malloc(max_size*sizeof(double));
    }
    
    printf("allocating space for rhs mat..\n");
    // make space for rhs vec
    rhs_mat = (double**)malloc(max_size*sizeof(double*));
    for(i=0; i<max_size; i++){ 
        rhs_mat[i] = (double*)malloc(max_size*sizeof(double));
    }

    // prefill initial solution with zeros and 
    // initialize rhs mat
    printf("initializing guess solution and rhs mats..\n");
    h = step_sizes[num_grid_levels-1];
    for(i=0; i<max_size; i++){
        for(j=0; j<max_size; j++){
            x = i*h;
            y = j*h; 
            if(i==0 || i==(max_size-1) || j==0 || j==(max_size-1))
                U_mat[i][j] = bndry_func(x,y);
            else
                U_mat[i][j] = 0;
            //printf("U_mat[%d][%d] = %f\n", i, j, U_mat[i][j]);
            rhs_mat[i][j] = rhs_func(x,y);
        }
    }


    // call some v-cycles    
    printf("calling mg_vcycle\n");
    for(i=0; i<num_v_cycles; i++)
        mg_vcycle(U_mat, rhs_mat, num_grid_levels-1, num_pre_iters, num_post_iters);
    printf("after mg_vcycle\n");

 
    // compute error norm and record solution in file
    FILE *fp;
    fp = fopen("data/output.txt","w");
    error_norm = 0;
    solution_norm = 0;
    percent_error = 0;
    true_norm = 0;
    max_size = grid_sizes[num_grid_levels-1];
    h = step_sizes[num_grid_levels-1];
    printf("printing output (max_size = %d) and (h = %f)\n", max_size, h);
    for(i=0; i<max_size; i++){
        for(j=0; j<max_size; j++){
            x = i*h;
            y = j*h; 
            true_norm = true_norm + sol_func(x,y)*sol_func(x,y);
            solution_norm = solution_norm + U_mat[i][j]*U_mat[i][j];
            error_norm = error_norm + pow(fabs( U_mat[i][j] - sol_func(x,y) ),2); 
            if(sol_func(x,y) > 1e-8){
                percent_error = percent_error + pow(fabs( U_mat[i][j] - sol_func(x,y) ),2)/(sol_func(x,y)*sol_func(x,y));
            }
            fprintf(fp, "%f  %f  %f  %f\n", x, y, U_mat[i][j], sol_func(x,y));
        }
    }
    true_norm = sqrt(true_norm);
    solution_norm = sqrt(solution_norm);
    error_norm = 100*sqrt(error_norm);
    percent_error = 100*percent_error;
    printf("\ntrue_norm = %f\n", true_norm);
    printf("solution_norm = %f\n", solution_norm);
    printf("\nerror_norm = %f %\n", error_norm); 
    printf("relative error norm = %f\n %", percent_error);
    fclose(fp);

    return 0;
}


// make space for arrays and set up grid level 
// data and solution and residual arrays
void startup(){
    int i, j, max_size;
    grid_sizes = (int *)malloc(num_grid_levels * sizeof(int));
    step_sizes = (double *)malloc(num_grid_levels * sizeof(double));
    for(i=0; i<num_grid_levels; i++){
        grid_sizes[i] = (int)pow(2,i) + 2;
        step_sizes[i] = 1.0/(grid_sizes[i] - 1);
        printf("grid_sizes[%d] = %d and step_sizes[%d] = %f\n", i, grid_sizes[i], i, step_sizes[i]); 
    }
}


void gauss_seidel_iter(double **U, double **rhs, int level){
    int i, j, n, iter_num, max_size;
    double x,y,h,temp_val;
    double **U_new;
    n = grid_sizes[level];
    h = step_sizes[level];
    printf("gauss_seidel iter with n = %d and h = %f--->\n", n, h);

    // set up new U
    max_size = n;
    U_new = (double**)malloc(max_size*sizeof(double*));
    for(i=0; i<max_size; i++){ 
        U_new[i] = (double*)malloc(max_size*sizeof(double));
    }

    // fill up new U with zeros
    for(i=0; i<max_size; i++){
        for(j=0; j<max_size; j++){
            U_new[i][j] = 0;
        }
    }

    // perform one cycle of gauss_seidel iteration
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            x = i*h;
            y = j*h;
            if(i==0 || i==(n-1) || j==0 || j==(n-1))
                U_new[i][j] = bndry_func(x,y);
            else{
                //U_new[i][j] = (0.25)*(U[i][j-1] + U[i-1][j] + U[i+1][j] + U[i][j+1] + (h*h)*rhs[i][j]);
                U_new[i][j] = (0.25)*(U_new[i][j-1] + U_new[i-1][j] + U[i+1][j] + U[i][j+1] + (h*h)*rhs[i][j]);
            }
        }
    }

    for(i=0; i<n; i++){
            for(j=0; j<n; j++){
                U[i][j] = U_new[i][j];
                //U_new[i][j] = 0;
            }
    }

    printf("end of gauss_seidel..\n");
}


void mg_vcycle(double **vec, double **rhs, int level, int num_pre_iters, int num_post_iters){


    int i,j,n, max_size;
    double h;
    double **rh;
    double **rH;
    double **eh;
    double **eH;

    n = grid_sizes[level];
    h = step_sizes[level];
    
    printf("in mg_vcycle with level = %d and n=%d and h=%f\n", level,n,h);

    if(level == 0){
        //for(i=0; i<num_pre_iters; i++)
            gauss_seidel_iter(vec,rhs,level);
    }   
    else if(level>0){

        // pre-smooth solution
        for(i=0; i<num_pre_iters; i++)
            gauss_seidel_iter(vec,rhs,level);


        //initialize residual mats (on fine and coarse grids)
        max_size = grid_sizes[level];
        rh = (double**)malloc(max_size*sizeof(double*));
        for(i=0; i<max_size; i++){ 
            rh[i] = (double*)malloc(max_size*sizeof(double));
        }
        for(i=0; i<max_size; i++){
            for(j=0; j<max_size; j++){
                rh[i][j] = 0;
            }
        }
        
        max_size = grid_sizes[level-1];
        rH = (double**)malloc(max_size*sizeof(double*));
        for(i=0; i<max_size; i++){ 
            rH[i] = (double*)malloc(max_size*sizeof(double));
        }
        for(i=0; i<max_size; i++){
            for(j=0; j<max_size; j++){
                rH[i][j] = 0;
            }
        }

        
        //compute residual
        printf("computing residual..\n");
        for(i=0; i<n; i++){
            for(j=0; j<n; j++){
                if( i==0 || i==(n-1) || j==0 || j==(n-1) )
                    rh[i][j] = 0;
                else
                    rh[i][j] = ((double)1.0/(h*h))*(vec[i+1][j] + vec[i-1][j] + vec[i][j+1] +
vec[i][j-1] - 4.0*vec[i][j]) + rhs[i][j];
            }
        }         
        printf("done computing residual..\n");
        
        // transfer residual to coarser grid
        fineToCoarse(rh, rH, level); 

        printf("form matrices eH and eh..\n");

        // form matrices eH and eh and fill them up with zeros
        max_size = grid_sizes[level-1]; 
        eH = (double**)malloc(max_size*sizeof(double*));
        for(i=0; i<max_size; i++){ 
            eH[i] = (double*)malloc(max_size*sizeof(double));
        }
        for(i=0; i<max_size; i++){
            for(j=0; j<max_size; j++){
                eH[i][j] = 0;
            }
        }

        max_size = grid_sizes[level]; 
        eh = (double**)malloc(max_size*sizeof(double*));
        for(i=0; i<max_size; i++){ 
            eh[i] = (double*)malloc(max_size*sizeof(double));
        }
        for(i=0; i<max_size; i++){
            for(j=0; j<max_size; j++){
                eh[i][j] = 0;
            }
        }


        // recursively call mg_vcycle
        mg_vcycle(eH, rH, level-1, num_pre_iters, num_post_iters);


        // transfer coarse error to fine grid
        coarseToFine(eH, eh, level); 

        // add correction to solution
        for(i=0; i<max_size; i++){
            for(j=0; j<max_size; j++){
                //printf("eh[i][j] = %f\n", eh[i][j]);
                vec[i][j] = vec[i][j] + eh[i][j];
            }
        }

        // post-smooth solution
        for(i=0; i<num_post_iters; i++)
            gauss_seidel_iter(vec,rhs,level);
    }

    
}


void coarseToFine(double **coarse_vec, double **fine_vec, int level){
    printf("in coarseToFine..\n");
    int i, j, i_f, j_f, num_coarse_elems, num_fine_elems, n, nc, nf;
    nc = grid_sizes[level-1];
    nf = grid_sizes[level];

    for(i=0; i<nf; i++){
        for(j=0; j<nf; j++){
            fine_vec[i][j] = 0;
        }
    }
   
    // step 1 
    for(i=1; i<nc; i++){
        for(j=1; j<nc; j++){
            i_f = 2*(i-1);
            j_f = 2*(j-1);
            fine_vec[i_f][j_f] = coarse_vec[i][j];
        }
    }

    // step 2 
    for(i=1; i<nc; i++){
        for(j=1; j<nc; j++){
            i_f = 2*(i-1);
            j_f = 2*(j-1);
            fine_vec[i_f][j_f] = 0.5*(coarse_vec[i][j] + fine_vec[i_f+1][j_f]);
        }
    }


    // step 3 
    for(i=1; i<nc; i++){
        for(j=1; j<nc; j++){
            i_f = 2*(i-1);
            j_f = 2*(j-1);
            fine_vec[i_f][j_f] = 0.5*(coarse_vec[i][j] + fine_vec[i_f][j_f+1]);
        }
    }


    // step 4 
    for(i=1; i<nc; i++){
        for(j=1; j<nc; j++){
            i_f = 2*(i-1);
            j_f = 2*(j-1);
            fine_vec[i_f][j_f] = 0.25*(coarse_vec[i][j] + fine_vec[i_f+1][j_f] +
fine_vec[i_f][j_f+1] + fine_vec[i_f+1][j_f+1]);
        }
    }

    printf("end of coarseToFine..\n");
}


void fineToCoarse(double **fine_vec, double **coarse_vec, int level){
    printf("in fineToCoarse..\n");
    int i, j, i_f, j_f, n, n1, n2, nc;
    n = grid_sizes[level];
    nc = grid_sizes[level-1];

    for(i=0; i<nc; i++){
        for(j=0; j<nc; j++){
            coarse_vec[i][j] = 0;
        }
    }

    for(i=1; i<nc; i++){
        for(j=1; j<nc; j++){
            i_f = 2*(i-1);
            j_f = 2*(j-1);
            if(i_f > 1 && j_f > 1){
                coarse_vec[i][j] = (1/16)*( 4*fine_vec[i_f][j_f] + 2*fine_vec[i_f+1][j_f] + 2*fine_vec[i_f-1][j_f] + 2*fine_vec[i_f][j_f+1] + 2*fine_vec[i_f][j_f-1] + fine_vec[i_f+1][j_f+1] + fine_vec[i_f+1][j_f-1] + fine_vec[i_f-1][j_f+1] + fine_vec[i_f-1][j_f-1] );
            }
            else{
                coarse_vec[i][j] = (1/4)*( fine_vec[i_f][j_f] + fine_vec[i_f+1][j_f] +
fine_vec[i_f][j_f+1] + fine_vec[i_f+1][j_f+1] );
            }
        }
    }

    printf("end of fineToCoarse..\n");
}


/* boundary condition */
double bndry_func(double x, double y){
    double res;
    if(problem_num == 1){
        res = 0;
    }
    else{
        res = (x*x + y*y);
    }
    return res;
}


/* rhs function */
double rhs_func(double x, double y){
    double res;
    double pi = 4.0*atan(1);
    if(problem_num == 1){
        res = 2*pi*pi*sin(pi*x)*sin(pi*y);  
    }
    else{
        res = -4.0;
    }
    return res;
}


/* analytical solution for comparison */
double sol_func(double x, double y){
    double res;
    double pi = 4.0*atan(1);
    if(problem_num == 1){
        res = sin(pi*x)*sin(pi*y);
    }
    else{
        res = x*x + y*y;
    }
    return res;
}

