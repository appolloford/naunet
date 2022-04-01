#include <stdio.h>
#include <stdexcept>

#include "naunet.h"
#include "naunet_data.h"
#include "naunet_macros.h"
#include "naunet_ode.h"
#include "naunet_timer.h"

int main() {
    int nsystem      = 2048;
    double spy       = 86400.0 * 365.0;

    double nH       = 2e4;
    double Tgas     = 15.0;

    NaunetData *data = new NaunetData[nsystem];
    for (int isys = 0; isys < nsystem; isys++) {
        data[isys].nH       = nH;
        data[isys].Tgas     = Tgas;
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
        // Set your initial abundance here
        throw std::runtime_error("Abundance has not been assigned");
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
