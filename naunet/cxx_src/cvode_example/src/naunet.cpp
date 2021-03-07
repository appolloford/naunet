#include "naunet.h"

// check_flag function is from the cvDiurnals_ky.c example from the CVODE package.
// Check function return value...
//   opt == 0 means SUNDIALS function allocates memory so check if
//            returned NULL pointer
//   opt == 1 means SUNDIALS function returns a flag so check if
//            flag >= 0
//   opt == 2 means function allocates memory so check if returned
//            NULL pointer
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL)
    {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return 1;
    }

    /* Check if flag < 0 */
    else if (opt == 1)
    {
        errflag = (int *)flagvalue;
        if (*errflag < 0)
        {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                    funcname, *errflag);
            return 1;
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL)
    {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return 1;
    }

    return 0;
}


Naunet::Naunet():
    m_atol(1e-20),
    m_rtol(1e-5)
{
    m_y = N_VNew_Serial((sunindextype)NSPECIES);
    m_cvode_mem = CVodeCreate(CV_BDF);
    m_ls = NULL;
};

Naunet::~Naunet()
{
    // N_VDestroy(m_y);
    CVodeFree(&m_cvode_mem);
    SUNLinSolFree(m_ls);
    // delete m_data;
}; 

int Naunet::initSolver()
{
    int flag;
    flag = CVodeInit(m_cvode_mem, fex, 0.0, m_y);
    if ( flag != CV_SUCCESS ) 
    {
        check_flag(&flag, "CVodeInit", 1);
    }
    CVodeSStolerances(m_cvode_mem, m_rtol, m_atol);
    if ( flag != CV_SUCCESS ) 
    {
        check_flag(&flag, "CVodeSStolerances", 1);
    }
    return 1;
};

int Naunet::solve(realtype *ab, realtype dt, UserData *data)
{
    int flag;
    N_VSetArrayPointer(ab, m_y);
    flag = CVodeReInit(m_cvode_mem, 0.0, m_y);
    if ( check_flag(&flag, "CVodeInit", 1) )
        return 1;
    flag = CVodeSetUserData(m_cvode_mem, data);
    if (check_flag(&flag, "CVodeSetUserData", 1))
        return 1;
    m_ls = SUNLinSol_SPGMR(m_y, 0, 0);
    if (check_flag((void *)m_ls, "SUNLinSol_SPGMR", 0))
        return 1;
    flag = CVSpilsSetLinearSolver(m_cvode_mem, m_ls);
    if (check_flag(&flag, "CVSpilsSetLinearSolver", 1))
        return 1;
    flag = CVSpilsSetJacTimes(m_cvode_mem, NULL, jtv);
    if (check_flag(&flag, "CVSpilsSetJacTimes", 1))
        return 1;

    realtype t0 = 0.0;
    flag = CVode(m_cvode_mem, dt, m_y, &t0, CV_NORMAL);
    for (int i=0; i<NSPECIES; i++) {
        ab[i] = NV_Ith_S(m_y, i);
    }
    return 0;
};



