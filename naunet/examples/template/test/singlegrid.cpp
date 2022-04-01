#include <stdio.h>
#include <stdexcept>

#include "naunet.h"
#include "naunet_data.h"
#include "naunet_macros.h"
#include "naunet_ode.h"
#include "naunet_timer.h"

int main() {
    double spy     = 86400.0 * 365.0;

    double nH       = 2e4;
    double Tgas     = 15.0;

    NaunetData data;
    data.nH       = nH;
    data.Tgas     = Tgas;

    Naunet naunet;
    naunet.Init();

#ifdef USE_CUDA
    naunet.Reset(1);
#endif

    double y[NEQUATIONS] = {0.0};
    for (int i = 0; i < NEQUATIONS; i++)
    {
        y[i] = 1.e-40;
    }
    // Set your initial abundance here
    throw std::runtime_error("Abundance has not been assigned");

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
