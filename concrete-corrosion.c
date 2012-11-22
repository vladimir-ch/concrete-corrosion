#include <stdlib.h>  /* standard library functions */
#include <stdio.h>   /* standard input/output functions */
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

/* CVODE header files */
#include <cvode/cvode.h>
#include <cvode/cvode_spgmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>

/* parameters */
#define BIOT_NUMBER       864.0
#define H                 2.5
#define C_BAR             1.0
#define PARTIAL_ORDER_P   1.0
#define PARTIAL_ORDER_Q   1.0

/* diffusivities */
#define D_1         864.0
#define D_2         0.0864
#define D_3         0.0864

/* rate constants */
#define ALPHA       7.2
#define BETA        0.84
#define K           1.0

/* boundary and initial conditions */
#define U1_D        0.1
#define U1_INIT     U1_D
#define U2_INIT     0.0
#define U3_INIT     0.00001
#define U4_INIT     0.0

/* geometry */
#define L_X         500.0
#define L_Y         10.0
#define N_X         256
#define N_Y         32
#define DX          (L_X/N_X)
#define DY          (L_Y/N_Y)
#define NEQ_X       (N_X+1)
#define NEQ_Y       (N_Y+1)
#define NEQ         (2*NEQ_X+2*NEQ_X*NEQ_Y-1)

/* time integration and output*/
#define T           1000.0
#define ABS_TOL     1.0e-9
#define REL_TOL     1.0e-9
#define N_OUT       100
#define DT          (T/N_OUT)

/* macros to access solution vector */
#define U1(u,xi)    u[(xi)*(2*NEQ_Y+2)-1]
#define U2(u,xi,yi) u[(xi)*(2*NEQ_Y+2)+2*(yi)]
#define U3(u,xi,yi) u[(xi)*(2*NEQ_Y+2)+2*(yi)+1]
#define U4(u,xi)    u[(xi)*(2*NEQ_Y+2)+2*NEQ_Y]

void   print_summary();
void   set_initial_profiles(N_Vector u_vec);
void*  initialize_cvode(N_Vector u_0_vec, CVRhsFn rhs_fn);
int    two_scale_corrosion_rhs(double t, N_Vector u_vec,
            N_Vector udot_vec, void* user_data);
double eta(double a, double b);
void   save_output(N_Vector u_vec);
void   print_final_cvode_stats(void* cvode_mem);

int main(void) {
    void*        cvode_mem;
    N_Vector     u_vec;
    double       t, tout;
    int          iout;

    u_vec = N_VNew_Serial((long)NEQ);
    if (u_vec == NULL) {
        exit(EXIT_FAILURE);
    }
    set_initial_profiles(u_vec);
    save_output(u_vec);

    cvode_mem = initialize_cvode(u_vec, two_scale_corrosion_rhs);
    
    print_summary();

    for (iout = 1; iout <= N_OUT; ++iout) {
        tout = iout*DT;
        if (CVode(cvode_mem, tout, u_vec, &t, CV_NORMAL) != 0) {
            exit(EXIT_FAILURE);
        }
        save_output(u_vec);
        printf("Step %5u / %5u\r", iout, N_OUT);
        fflush(stdout);
    }
    print_final_cvode_stats(cvode_mem);
    N_VDestroy_Serial(u_vec);
    CVodeFree(&cvode_mem);
    return EXIT_SUCCESS;
}

void set_initial_profiles(N_Vector u_vec) {
    double* u = NV_DATA_S(u_vec);
    int xi, yi;

    for (xi = 0; xi <= N_X; ++xi) {
        if (xi != 0) {
            U1(u,xi) = U1_INIT;
        }
        for (yi = 0; yi <= N_Y; ++yi) {
            U2(u,xi,yi) = U2_INIT;
            U3(u,xi,yi) = U3_INIT;
        }
        U4(u,xi) = U4_INIT;
    }
}

void* initialize_cvode(N_Vector u_0_vec, CVRhsFn rhs_fn) {
    int flag;
    
    void* cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (cvode_mem == NULL) { exit(EXIT_FAILURE); }
    flag = CVodeInit(cvode_mem, rhs_fn, 0.0, u_0_vec);
    if (flag != 0 ) { exit(EXIT_FAILURE); }
    flag = CVodeSStolerances(cvode_mem, REL_TOL, ABS_TOL);
    if (flag != 0 ) { exit(EXIT_FAILURE); }
    flag = CVodeSetMaxNumSteps(cvode_mem, 1000000L); 
    if (flag != 0 ) { exit(EXIT_FAILURE); }
    flag = CVSpgmr(cvode_mem, PREC_NONE, 30);
    if (flag != 0 ) { exit(EXIT_FAILURE); }
    return cvode_mem;
}

int two_scale_corrosion_rhs(double t, N_Vector u_vec,
        N_Vector udot_vec, void *user_data)
{
    double* u = NV_DATA_S(u_vec);
    double* udot = NV_DATA_S(udot_vec);
    double  u1, u2, u3;
    double  BHu1u2, etau3u4;
    double  u1_left, u1_right;
    double  u2_left, u2_right;
    double  u3_left, u3_right;
    int     xi, yi;

    for (xi = 0; xi <= N_X; ++xi) {
        if (xi != 0) {
            u1 = U1(u,xi);
        } else {
            u1 = U1_D;
        }
        BHu1u2 = BIOT_NUMBER*(H*u1 - U2(u,xi,0));
        if (xi != 0) {
            if (xi != 1) {
                u1_left = U1(u,xi-1);
            } else {
                u1_left = U1_D;
            }
            u1_right = (xi != N_X) ? U1(u,xi+1) : u1_left;
            U1(udot,xi) = D_1*(u1_left - 2.0*u1 + u1_right)/(DX*DX)
              - BHu1u2;
        }
        for (yi = 0; yi <= N_Y; ++yi) {
            if (yi != 0) {
                u2_left = U2(u,xi,yi-1);
                u3_left = U3(u,xi,yi-1);
            } else {
                u2_left = U2(u,xi,1) + 2.0*DY/D_2*BHu1u2;
                u3_left = U3(u,xi,1);
            }
            u2 = U2(u,xi,yi);
            u3 = U3(u,xi,yi);
            if (yi != N_Y) {
                u2_right = U2(u,xi,yi+1);
                u3_right = U3(u,xi,yi+1);
            } else {
                etau3u4 = eta(u3, U4(u,xi));
                u2_right = u2_left;
                u3_right = u3_left - 2.0*DY/D_3*etau3u4;
            }
            U2(udot,xi,yi) = D_2*(u2_left - 2.0*u2 + u2_right)/(DY*DY)
              - ALPHA*u2 + BETA*u3;
            U3(udot,xi,yi) = D_3*(u3_left - 2.0*u3 + u3_right)/(DY*DY)
              + ALPHA*u2 - BETA*u3;
        }
        U4(udot,xi) = etau3u4;
    }

    return 0;
}

double eta(double a, double b) {
    if (a > 0.0 && b < C_BAR) {
        return K*pow(a, PARTIAL_ORDER_P)*pow(C_BAR - b, PARTIAL_ORDER_Q);
    } else {
        return 0.0;
    }
}

void print_summary() {
    printf("Omega = (0,%f)\n", L_X);
    printf("Y = (0,%f)\n", L_Y);
    printf("grid size = %u x %u\n", NEQ_X, NEQ_Y);
    printf("step sizes: h_x = %f, h_y = %f\n", DX, DY);
    printf("d_1 = %f, d_2 = %f, d_3 = %f\n", D_1, D_2, D_3);
    printf("alpha = %f, beta = %f, k = %f\n", ALPHA, BETA, K);
    printf("Biot = %f, c_bar = %f, H = %f, p = %f, q = %f\n",
        BIOT_NUMBER, C_BAR, H, PARTIAL_ORDER_P, PARTIAL_ORDER_Q);
    printf("u1_D = %f\n", U1_D);
    printf("u1_0 = %f, u2_0 = %f, u3_0 = %f, u4_0 = %f\n",
        U1_INIT, U2_INIT, U3_INIT, U4_INIT);
    printf("time interval = (0,%f)\n", T);
    printf("relative tolerance = %g\nabsolute tolerance = %g\n\n",
        REL_TOL, ABS_TOL);
}

void save_output(N_Vector u_vec) {
    static int  count = 0;
    static char output_dir[256];
    char        filename[512];
    int         xi, yi;
    FILE*       fout;
    double*     u = NV_DATA_S(u_vec);
    
    if (count == 0) {
        time_t     rawtime;
        struct tm* timeinfo;
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(output_dir, (size_t)255, "result_%Y%m%d_%H%m%S", timeinfo);
        if (0 != mkdir(output_dir, 0755)) {
            perror("Failed to create output directory");
            exit(EXIT_FAILURE);
        } else {
            printf("Writing results to directory: %s\n", output_dir);
        }
    }

    sprintf(filename, "%s/u1_%05u.dat", output_dir, count); 
    fout = fopen(filename, "wt");
    if (fout == NULL) {
        perror("Failed to create output file for u1");
        exit(EXIT_FAILURE);
    }
    fprintf(fout, "%f %f\n", 0.0, U1_D);
    for (xi = 1; xi <= N_X; ++xi) {
        fprintf(fout, "%f %f\n", xi*DX, U1(u,xi));
    }
    fclose(fout);

    sprintf(filename, "%s/u2_%05u.dat", output_dir, count); 
    fout = fopen(filename, "wt");
    if (fout == NULL) {
        perror("Failed to create output file for u2");
        exit(EXIT_FAILURE);
    }
    for (xi = 0; xi <= N_X; ++xi) {
        for (yi = 0; yi <= N_Y; ++yi) {
            fprintf(fout, "%f %f %f\n", xi*DX, yi*DY, U2(u,xi,yi));
        }
        fprintf(fout, "\n");
    }
    fclose(fout);

    sprintf(filename, "%s/u3_%05u.dat", output_dir, count); 
    fout = fopen(filename, "wt");
    if (fout == NULL) {
        perror("Failed to create output file for u3");
        exit(EXIT_FAILURE);
    }
    for (xi = 0; xi <= N_X; ++xi) {
        for (yi = 0; yi <= N_Y; ++yi) {
            fprintf(fout, "%f %f %f\n", xi*DX, yi*DY, U3(u,xi,yi));
        }
        fprintf(fout, "\n");
    }
    fclose(fout);

    sprintf(filename, "%s/u4_%05u.dat", output_dir, count); 
    fout = fopen(filename, "wt");
    if (fout == NULL) {
        perror("Failed to create output file for u4");
        exit(EXIT_FAILURE);
    }
    for (xi = 0; xi <= N_X; ++xi) {
        fprintf(fout, "%f %f\n", xi*DX, U4(u, xi));
    }
    fclose(fout);

    ++count;
}

void print_final_cvode_stats(void* cvode_mem) {
    long int lenrw, leniw ;
    long int lenrwLS, leniwLS;
    long int nst, nfe, nsetups, nni, ncfn, netf;
    long int nli, npe, nps, ncfl, nfeLS;

    /* CVODE solver stats */
    CVodeGetWorkSpace(cvode_mem, &lenrw, &leniw);
    CVodeGetNumSteps(cvode_mem, &nst);
    CVodeGetNumRhsEvals(cvode_mem, &nfe);
    CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    CVodeGetNumErrTestFails(cvode_mem, &netf);
    CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);

    /* Spils (linear solver) stats */
    CVSpilsGetWorkSpace(cvode_mem, &lenrwLS, &leniwLS);
    CVSpilsGetNumLinIters(cvode_mem, &nli);
    CVSpilsGetNumPrecEvals(cvode_mem, &npe);
    CVSpilsGetNumPrecSolves(cvode_mem, &nps);
    CVSpilsGetNumConvFails(cvode_mem, &ncfl);
    CVSpilsGetNumRhsEvals(cvode_mem, &nfeLS);

    printf("\n\nFinal CVODE statistics:\n\n");
    printf("Real workspace size = %5ld\n", lenrw);
    printf("Integer workspace size = %5ld\n", leniw);
    printf("Linear solver real workspace size = %5ld\n", lenrwLS);
    printf("Linear solver integer workspace size = %5ld\n", leniwLS);
    printf("Number of steps taken = %5ld\n", nst);
    printf("Number of RHS evaluations = %5ld\n", nfe);
    printf("Number of RHS evaluations in linear solver = %5ld\n", nfeLS);
    printf("Number of nonlinear solver iterations = %5ld\n", nni);
    printf("Number of linear solver iterations = %5ld\n", nli);
    printf("Number of linear solver setups = %5ld\n", nsetups);
    printf("Number of error test fails = %5ld\n", netf);
    printf("Number of preconditioner evaluations = %5ld\n", npe);
    printf("Number of preconditioner solves = %5ld\n", nps);
    printf("Number of nonlinear convergence fails = %5ld\n", ncfn);
    printf("Number of linear solver convergence fails = %5ld\n\n", ncfl);
}

