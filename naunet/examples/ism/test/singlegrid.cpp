#include <stdio.h>

#include "naunet.h"
#include "naunet_userdata.h"
#include "naunet_macros.h"
#include "naunet_ode.h"
#include "naunet_timer.h"

int main()
{

    double spy = 86400.0 * 365.0;

    double nH = 1e4;
    double zeta_cr = 1e-16;
    double zeta_xr = 0.0;
    double Tgas = 15.0;
    double Tdust = 15.0;
    double Av = 100.0;
    double G0 = 4.0;
    double rG = 1e-5;
    double omega = 0.5;
    double barr = 1.5e-8;
    double sites = 1e15;
    double hop = 0.3;
    double nMono = 2.0;
    double duty = 3.16e-19;
    double Tcr = 70.0;
    double branch = 1e-2;

    UserData data;
    data.nH = 1e4;
    data.zeta_cr = 1e-16;
    data.zeta_xr = 0.0;
    data.Tgas = 15.0;
    data.Tdust = 15.0;
    data.Av = 100.0;
    data.G0 = 4.0;
    data.rG = 1e-5;
    data.omega = 0.5;
    data.barr = 1.5e-8;
    data.sites = 1e15;
    data.hop = 0.3;
    data.nMono = 2.0;
    data.duty = 3.16e-19;
    data.Tcr = 70.0;
    data.branch = 1e-2;

    Naunet naunet;
    naunet.initSolver();

#ifdef USE_CUDA
    naunet.resetSolver(1);
#endif

    double y[NSPECIES];
    for (int i = 0; i < NSPECIES; i++)
    {
        y[i] = 1.e-40;
    }
    y[IDX_H2I] = 0.5 * nH;
    y[IDX_HI] = 5.0e-5 * nH;
    y[IDX_HeI] = 9.75e-2 * nH;
    y[IDX_NI] = 7.5e-5 * nH;
    y[IDX_OI] = 3.2e-4 * nH;
    y[IDX_CI] = 1.4e-4 * nH;
    y[IDX_SI] = 8.0e-8 * nH;
    y[IDX_SiI] = 8.0e-9 * nH;
    y[IDX_NaI] = 2.0e-9 * nH;
    y[IDX_MgI] = 7.0e-9 * nH;
    y[IDX_FeI] = 3.0e-9 * nH;
    y[IDX_ClI] = 4.0e-9 * nH;
    y[IDX_FI] = 2.0e-8 * nH;
    y[IDX_GRAIN0I] = 1.3e-12 * nH;

    FILE *fbin = fopen("evolution_singlegrid.bin", "w");
    FILE *ftxt = fopen("evolution_singlegrid.txt", "w");
    FILE *ttxt = fopen("time_singlegrid.txt", "w");
#ifdef NAUNET_DEBUG
    printf("Initialization is done. Start to evolve.\n");
    FILE *rtxt = fopen("reactionrates.txt", "w");
    double rates[NREACTIONS];
#endif

    double logtstart = 3.0, logtend = 8.0;
    double dtyr = 0.0, time = 0.0;
    for (double logtime = logtstart; logtime < logtend; logtime += 0.1)
    {

#ifdef NAUNET_DEBUG
        calculate_rates(rates, y, &data);
        for (int j = 0; j < NREACTIONS; j++)
        {
            fprintf(rtxt, "%13.7e ", rates[j]);
        }
        fprintf(rtxt, "\n");
#endif

        dtyr = pow(10.0, logtime) - time;
        time += dtyr;

        fwrite(&time, sizeof(double), 1, fbin);
        fwrite(y, sizeof(double), NSPECIES, fbin);

        fprintf(ftxt, "%13.7e ", time);
        for (int j = 0; j < NSPECIES; j++)
        {
            fprintf(ftxt, "%13.7e ", y[j]);
        }
        fprintf(ftxt, "\n");

        Timer timer;
        timer.start();
        naunet.solve(y, dtyr * spy, &data);
        timer.stop();
        float duration = (float)timer.elapsed() / 1e6;
        fprintf(ttxt, "%8.5e \n", duration);
        // printf("Time = %13.7e yr, elapsed: %8.5e sec\n", time[i + 1], duration);
    }

    fclose(fbin);
    fclose(ftxt);
    fclose(ttxt);
#ifdef NAUNET_DEBUG
    fclose(rtxt);
#endif

    return 0;
}