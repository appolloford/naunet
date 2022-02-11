#include <stdio.h>

#include "naunet.h"
#include "naunet_data.h"
#include "naunet_macros.h"
#include "naunet_ode.h"
#include "naunet_timer.h"

int main() {
    int nsystem      = 2048;
    double spy       = 86400.0 * 365.0;

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
    double eb_crd   = 1.21e3;
    double eb_h2d   = 1.21e3;
    double eb_uvd   = 1.00e4;
    double crdeseff = 1.0e5;
    double h2deseff = 1.0e-2;
    double uvcreff  = 1.0e-3;

    NaunetData *data = new NaunetData[nsystem];
    for (int isys = 0; isys < nsystem; isys++) {
        data[isys].nH       = nH;
        data[isys].zeta     = zeta;
        data[isys].Tgas     = Tgas;
        data[isys].Av       = Av;
        data[isys].omega    = omega;
        data[isys].G0       = G0;
        data[isys].rG       = rG;
        data[isys].gdens    = gdens;
        data[isys].sites    = sites;
        data[isys].fr       = fr;
        data[isys].opt_thd  = opt_thd;
        data[isys].opt_crd  = opt_crd;
        data[isys].opt_uvd  = opt_uvd;
        data[isys].opt_h2d  = opt_h2d;
        data[isys].eb_crd   = eb_crd;
        data[isys].eb_h2d   = eb_h2d;
        data[isys].eb_uvd   = eb_uvd;
        data[isys].crdeseff = crdeseff;
        data[isys].h2deseff = h2deseff;
        data[isys].uvcreff  = uvcreff;
    }

    Naunet naunet;
    if (naunet.Init() == NAUNET_FAIL) {
        printf("Initialize Fail\n");
    }

    if (naunet.Reset(nsystem) == NAUNET_FAIL) {
        printf("Reset Fail\n");
    }

    double *y = new double[nsystem * NEQUATIONS];
    for (int isys = 0; isys < nsystem; isys++) {
        for (int i = 0; i < NEQUATIONS; i++) {
            y[isys * NEQUATIONS + i] = 1.e-40;
        }
        y[isys * NEQUATIONS + IDX_H2I] = 0.5 * nH;
        y[isys * NEQUATIONS + IDX_HI]  = 5.0e-5 * nH;
        y[isys * NEQUATIONS + IDX_HeI] = 9.75e-2 * nH;
        y[isys * NEQUATIONS + IDX_NI]  = 7.5e-5 * nH;
        y[isys * NEQUATIONS + IDX_OI]  = 1.8e-4 * nH;
        y[isys * NEQUATIONS + IDX_COI] = 1.4e-4 * nH;
        y[isys * NEQUATIONS + IDX_SI]  = 8.0e-8 * nH;
        y[isys * NEQUATIONS + IDX_SiI] = 8.0e-9 * nH;
        y[isys * NEQUATIONS + IDX_MgI] = 7.0e-9 * nH;
        y[isys * NEQUATIONS + IDX_ClI] = 4.0e-9 * nH;
    }

    FILE *fbin = fopen("evolution_multiplegrid.bin", "w");
    FILE *ftxt = fopen("evolution_multiplegrid.txt", "w");
    FILE *ttxt = fopen("time_parallel.txt", "w");

#ifdef NAUNET_DEBUG
    printf("Initialization is done. Start to evolve.\n");
    // FILE *rtxt = fopen("reactionrates.txt", "w");
    // double rates[NREACTIONS];
#endif

    double logtstart = 2.0, logtend = 4.0;
    double dtyr = 0.0, time = 0.0;
    for (double logtime = logtstart; logtime < logtend; logtime += 0.1) {
#ifdef NAUNET_DEBUG
        // EvalRates only receive one system as input, disabled in parallel test
        // EvalRates(rates, y, data);
        // for (int j = 0; j < NREACTIONS; j++) {
        //     fprintf(rtxt, "%13.7e ", rates[j]);
        // }
        // fprintf(rtxt, "\n");
#endif

        dtyr = pow(10.0, logtime) - time;

        for (int isys = 0; isys < nsystem; isys++) {
            fwrite((double *)&isys, sizeof(double), 1, fbin);
            fwrite(&time, sizeof(double), 1, fbin);
            fwrite(&y[isys * NEQUATIONS], sizeof(double), NEQUATIONS, fbin);

            fprintf(ftxt, "%13.7e ", (double)isys);
            fprintf(ftxt, "%13.7e ", time);
            for (int j = 0; j < NEQUATIONS; j++) {
                fprintf(ftxt, "%13.7e ", y[isys * NEQUATIONS + j]);
            }
            fprintf(ftxt, "\n");
        }

        Timer timer;
        timer.start();
        naunet.Solve(y, dtyr * spy, data);
        timer.stop();

        time += dtyr;

        // float duration = (float)timer.elapsed() / 1e6;
        double duration = timer.elapsed();
        fprintf(ttxt, "%8.5e \n", duration);
        printf("Time = %13.7e yr, elapsed: %8.5e sec\n", time, duration);
    }

    for (int isys = 0; isys < nsystem; isys++) {
        fwrite((double *)&isys, sizeof(double), 1, fbin);
        fwrite(&time, sizeof(double), 1, fbin);
        fwrite(&y[isys * NEQUATIONS], sizeof(double), NEQUATIONS, fbin);

        fprintf(ftxt, "%13.7e ", (double)isys);
        fprintf(ftxt, "%13.7e ", time);
        for (int j = 0; j < NEQUATIONS; j++) {
            fprintf(ftxt, "%13.7e ", y[isys * NEQUATIONS + j]);
        }
        fprintf(ftxt, "\n");
    }

    fclose(fbin);
    fclose(ftxt);
    fclose(ttxt);

#ifdef NAUNET_DEBUG
    // fclose(rtxt);
#endif

    if (naunet.Finalize() == NAUNET_FAIL) {
        printf("Finalize Fail\n");
    }

    delete[] data;
    delete[] y;

    return 0;
}
