#ifndef __NAUNET_H__
#define __NAUNET_H__

#include <cvode/cvode.h>                // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>     // access to serial N_Vector
#include <sunlinsol/sunlinsol_spgmr.h>  // access to SPGMR SUNLinearSolver
#include <cvode/cvode_spils.h>          // access to CVSpils interface
#include <sunmatrix/sunmatrix_dense.h>  // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h>  // access to dense SUNLinearSolver
#include <sundials/sundials_dense.h>    // use generic dense solver in precond
#include <sundials/sundials_types.h>    // defs. of realtype, sunindextype
#include <sundials/sundials_math.h>     // contains the macros ABS, SUNSQR, EXP
#include <sunmatrix/sunmatrix_sparse.h> // access to sparse SUNMatrix
#include <sunlinsol/sunlinsol_klu.h>    // access to KLU sparse direct solver

#include "naunet_userdata.h"
#include "naunet_constants.h"
#include "naunet_ode.h"

class Naunet
{

private:
    // UserData *m_data;
    N_Vector m_y;
    SUNMatrix m_a;

    realtype m_atol;
    realtype m_rtol;
    void *m_cvode_mem;
    SUNLinearSolver m_ls;

public:
    Naunet();
    ~Naunet();
    int initSolver();
    int solve(realtype *ab, realtype dt, UserData *data);
    // Naunet(){
    //     // m_data = new UserData();
    //     m_y = N_VNew_Serial((sunindextype)NSPECIES);
    //     m_cvode_mem = CVodeCreate(CV_BDF);
    //     m_ls = NULL;
    // };

    // int initSolver(){
    //     int flag;
    //     flag = CVodeInit(m_cvode_mem, fex, 0.0, m_y);
    //     if ( check_flag(&flag, "CVodeInit", 1) )
    //         return 1;
    //     flag = CVodeSStolerances(m_cvode_mem, m_rtol, m_atol);
    //     if ( check_flag(&flag, "CVodeSStolerances", 1) )
    //         return 1;
    //     return 0;
    // };

    // int solve(realtype *ab, realtype dt, UserData *data){
    //     int flag;
    //     N_VSetArrayPointer(ab, m_y);
    //     flag = CVodeReInit(m_cvode_mem, 0.0, m_y);
    //     if ( check_flag(&flag, "CVodeInit", 1) )
    //         return 1;
    //     flag = CVodeSetUserData(m_cvode_mem, data);
    //     if (check_flag(&flag, "CVodeSetUserData", 1))
    //         return 1;
    //     m_ls = SUNLinSol_SPGMR(m_y, 0, 0);
    //     if (check_flag((void *)m_ls, "SUNLinSol_SPGMR", 0))
    //         return 1;
    //     flag = CVSpilsSetLinearSolver(m_cvode_mem, m_ls);
    //     if (check_flag(&flag, "CVSpilsSetLinearSolver", 1))
    //         return 1;
    //     flag = CVSpilsSetJacTimes(m_cvode_mem, NULL, jtv);
    //     if (check_flag(&flag, "CVSpilsSetJacTimes", 1))
    //         return 1;

    //     realtype t0 = 0.0;
    //     flag = CVode(m_cvode_mem, dt, m_y, &t0, CV_NORMAL);
    //     for (int i=0; i<NSPECIES; i++) {
    //         ab[i] = NV_Ith_S(m_y, i);
    //     }
    //     return 0;
    // };

    // ~Naunet(){
    //     // N_VDestroy(m_y);
    //     CVodeFree(&m_cvode_mem);
    //     SUNLinSolFree(m_ls);
    //     // delete m_data;
    // };
};

#endif
