#include <stdio.h>

#include "naunet.h"
#include "naunet_data.h"
#include "naunet_macros.h"
#include "naunet_ode.h"
#include "naunet_timer.h"

int main() {
    double spy     = 86400.0 * 365.0;

    double nH       = 2e4;
    double zeta     = 1.3e-17;
    double Tgas     = 15.0;
    double Av       = 10.0;
    double omega    = 0.5;
    double G0       = 1.0;
    double rG       = 1e-5;
    double gdens    = 7.6394373e-13 * nH;
    double sites    = 1.5e15;
    double fr       = 1.0;
    double opt_thd  = 1.0;
    double opt_crd  = 1.0;
    double opt_uvd  = 1.0;
    double opt_h2d  = 1.0;
    double crdeseff = 1.0e5;
    double h2deseff = 1.0e-2;
    double uvcreff  = 1.0e-3;

    NaunetData data;
    data.nH       = nH;
    data.zeta     = zeta;
    data.Tgas     = Tgas;
    data.Av       = Av;
    data.omega    = omega;
    data.G0       = G0;
    data.rG       = rG;
    data.gdens    = gdens;
    data.sites    = sites;
    data.fr       = fr;
    data.opt_thd  = opt_thd;
    data.opt_crd  = opt_crd;
    data.opt_uvd  = opt_uvd;
    data.opt_h2d  = opt_h2d;
    data.crdeseff = crdeseff;
    data.h2deseff = h2deseff;
    data.uvcreff  = uvcreff;

    Naunet naunet;
    naunet.Init();

#ifdef USE_CUDA
    naunet.Reset(1);
#endif

    double y[NEQUATIONS] = {0.0};
    // for (int i = 0; i < NEQUATIONS; i++)
    // {
    //     y[i] = 1.e-40;
    // }
    y[IDX_H2I]           = 0.5 * nH;
    y[IDX_HI]            = 5.0e-5 * nH;
    y[IDX_HeI]           = 9.75e-2 * nH;
    y[IDX_NI]            = 7.5e-5 * nH;
    y[IDX_OI]            = 1.8e-4 * nH;
    y[IDX_COI]           = 1.4e-4 * nH;
    y[IDX_SI]            = 8.0e-8 * nH;
    y[IDX_SiI]           = 8.0e-9 * nH;
    y[IDX_MgI]           = 7.0e-9 * nH;
    y[IDX_ClI]           = 4.0e-9 * nH;

    FILE *fbin           = fopen("evolution_singlegrid.bin", "w");
    FILE *ftxt           = fopen("evolution_singlegrid.txt", "w");
    FILE *ttxt           = fopen("time_singlegrid.txt", "w");
#ifdef NAUNET_DEBUG
    printf("Initialization is done. Start to evolve.\n");
    FILE *rtxt = fopen("reactionrates.txt", "w");
    double rates[NREACTIONS] = {0.0};
#endif

    double logtstart = 3.0, logtend = 7.0;
    double dtyr = 0.0, time = 0.0;
    for (double logtime = logtstart; logtime < logtend; logtime += 0.1) {
#ifdef NAUNET_DEBUG
        EvalRates(rates, y, &data);
        for (int j = 0; j < NREACTIONS; j++)
        {
            fprintf(rtxt, "%13.7e ", rates[j]);
        }
        fprintf(rtxt, "\n");
#endif

        dtyr = pow(10.0, logtime) - time;

        fwrite(&time, sizeof(double), 1, fbin);
        fwrite(y, sizeof(double), NEQUATIONS, fbin);

        fprintf(ftxt, "%13.7e ", time);
        for (int j = 0; j < NEQUATIONS; j++) {
            fprintf(ftxt, "%13.7e ", y[j]);
        }
        fprintf(ftxt, "\n");

        Timer timer;
        timer.start();
        naunet.Solve(y, dtyr * spy, &data);
        timer.stop();

        time += dtyr;

        // float duration = (float)timer.elapsed() / 1e6;
        double duration = timer.elapsed();
        fprintf(ttxt, "%8.5e \n", duration);
        printf("Time = %13.7e yr, elapsed: %8.5e sec\n", time, duration);
    }

    // save the final results
    fwrite(&time, sizeof(double), 1, fbin);
    fwrite(y, sizeof(double), NEQUATIONS, fbin);

    fprintf(ftxt, "%13.7e ", time);
    for (int j = 0; j < NEQUATIONS; j++) {
        fprintf(ftxt, "%13.7e ", y[j]);
    }
    fprintf(ftxt, "\n");

    fclose(fbin);
    fclose(ftxt);
    fclose(ttxt);
#ifdef NAUNET_DEBUG
    fclose(rtxt);
#endif

    naunet.Finalize();

    return 0;
}
