/*
A simple example using the CVODE library to solve a simple 2d ODE, treating it
as a stiff system. Added the use of user pointer which is a critical building
block for solving more complex problems.
*/

#include <iostream>
#include <vector>
#include <cvode/cvode.h>               // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>    // access to serial N_Vector
#include <sunlinsol/sunlinsol_spgmr.h> //access to SPGMR SUNLinearSolver
#include <cvode/cvode_spils.h>         // access to CVSpils interface
#include <sundials/sundials_dense.h>   // use generic dense solver in precond
#include <sundials/sundials_types.h>   // defs. of realtype, sunindextype
#include <sundials/sundials_math.h>    // contains the macros ABS, SUNSQR, EXP

#include "naunet_userdata.h"
#include "naunet_constants.h"
#include "fex.h"
#include "jtv.h"

static int check_flag(void *flagvalue, const char *funcname, int opt);

int main()
{
    // Setup User Data Pointer
    UserData *data = new UserData();
    data->Tgas = 15.0;

    int flag;               // For checking if functions have run properly
    realtype abstol = 1e-5; // real tolerance of system
    realtype reltol = 1e-5; // absolute tolerance of system

    // 1. Initialize parallel or multi-threaded environment, if appropriate.
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    // 2. Defining the length of the problem.
    // ---------------------------------------------------------------------------
    sunindextype N = NSPECIES;
    // ---------------------------------------------------------------------------

    // 3. Set vector of initial values.
    // ---------------------------------------------------------------------------
    N_Vector y; // Problem vector.
    y = N_VNew_Serial(N);
    NV_Ith_S(y, 0) = 0.4;
    NV_Ith_S(y, 1) = 0.4;
    NV_Ith_S(y, 2) = 0.1;
    NV_Ith_S(y, 3) = 0.1;
    // ---------------------------------------------------------------------------

    // 4. Create CVODE Object.
    // ---------------------------------------------------------------------------
    void *cvode_mem = NULL; // Problem dedicated memory.
    cvode_mem = CVodeCreate(CV_BDF);
    // ---------------------------------------------------------------------------

    // 5. Initialize CVODE solver.
    // ---------------------------------------------------------------------------
    realtype t0 = 0; // Initiale value of time.
    flag = CVodeInit(cvode_mem, fex, t0, y);
    if (check_flag(&flag, "CVodeSetUserData", 1))
        return (1);
    // ---------------------------------------------------------------------------

    // 6. Specify integration tolerances.
    // ---------------------------------------------------------------------------
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1))
        return (1);
    // ---------------------------------------------------------------------------

    // 7. Set Optional inputs.
    // ---------------------------------------------------------------------------
    /* Set the pointer to user-defined data */
    flag = CVodeSetUserData(cvode_mem, data);
    // ---------------------------------------------------------------------------

    // 8. Create Matrix Object.
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    // 9. Create Linear Solver Object.
    // ---------------------------------------------------------------------------
    SUNLinearSolver LS;
    // Here we chose one of the possible linear solver modules. SUNSPMR is an
    // iterative solver that is designed to be compatible with any nvector
    // implementation (serial, threaded, parallel,
    // user-supplied)that supports a minimal subset of operations.
    LS = SUNLinSol_SPGMR(y, 0, 0);
    if (check_flag((void *)LS, "SUNSPGMR", 0))
        return (1);
    // ---------------------------------------------------------------------------

    // 10. Set linear solver optional inputs.
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    // 11. Attach linear solver module.
    // ---------------------------------------------------------------------------
    // CVSpilsSetLinearSolver is for iterative linear solvers.
    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
    if (check_flag(&flag, "CVSpilsSetLinearSolver", 1))
        return 1;
    // ---------------------------------------------------------------------------

    // 12. Set linear solver interface optional inputs.
    // ---------------------------------------------------------------------------
    // Sets the jacobian-times-vector function.
    flag = CVSpilsSetJacTimes(cvode_mem, NULL, jtv);
    if (check_flag(&flag, "CVSpilsSetJacTimes", 1))
        return (1);

    // ---------------------------------------------------------------------------

    // 13. Specify rootfinding problem.
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    // 14. Advance solution in time.
    // ---------------------------------------------------------------------------
    // Have the solution advance over time, but stop to log 100 of the steps.
    int print_steps = 100;
    realtype tout;
    realtype end_time = 50;
    realtype step_length = 0.5;
    realtype t = 0;
    // loop over output points, call CVode, print results, test for error
    for (tout = step_length; tout <= end_time; tout += step_length)
    {
        flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        std::cout << "t: " << t;
        std::cout << "\ny:";
        N_VPrint_Serial(y);
        if (check_flag(&flag, "CVode", 1))
            break;
    }
    // ---------------------------------------------------------------------------

    // 15. Get optional outputs.
    // ---------------------------------------------------------------------------
    // ---------------------------------------------------------------------------

    // 16. Deallocate memory for solution vector.
    // ---------------------------------------------------------------------------
    N_VDestroy(y);
    // ---------------------------------------------------------------------------

    // 17. Free solver memory.
    // ---------------------------------------------------------------------------
    CVodeFree(&cvode_mem);
    // ---------------------------------------------------------------------------

    // 18. Free linear solver and matrix memory.
    // ---------------------------------------------------------------------------
    SUNLinSolFree(LS);
    delete data;
    // Remember to free the user data memory.
    // ---------------------------------------------------------------------------

    return 0;
}

// check_flag function is from the cvDiurnals_ky.c example from the CVODE
// package.
/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL)
    {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return (1);
    }

    /* Check if flag < 0 */
    else if (opt == 1)
    {
        errflag = (int *)flagvalue;
        if (*errflag < 0)
        {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                    funcname, *errflag);
            return (1);
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL)
    {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return (1);
    }

    return (0);
}
