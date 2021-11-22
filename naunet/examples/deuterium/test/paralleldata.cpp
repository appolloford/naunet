#include <stdio.h>

#include "naunet.h"
#include "naunet_data.h"
#include "naunet_macros.h"
#include "naunet_ode.h"
#include "naunet_timer.h"

int main() {
    int nsystem               = 4096;
    double spy                = 86400.0 * 365.0;
    double pi                 = 3.14159265;
    double rD                 = 1.0e-5;
    double rhoD               = 3.0;
    double DtoGM              = 7.09e-3;
    double amH                = 1.66043e-24;
    double nH                 = 1e5;
    double OPRH2              = 0.1;

    double rawdata[4096][133] = {0.0};

    FILE *fab                 = fopen("grids.dat", "r");

    for (int i = 0; i < 4096; i++) {
        if (feof(fab)) break;

        for (int j = 0; j < 133; j++) {
            fscanf(fab, "%lf", &(rawdata[i][j]));
        }
    }

    fclose(fab);

#ifdef NAUNET_DEBUG
    printf("First row of input data\n");
    for (int i = 0; i < 133; i++) {
        printf("%13.7e ", rawdata[0][i]);
    }
    printf("\n");
#endif

    NaunetData *data = new NaunetData[nsystem];
    for (int isys = 0; isys < nsystem; isys++) {
        data[isys].nH          = rawdata[isys][0];
        data[isys].Tgas        = rawdata[isys][1];
        // data[isys].Tgas = 15.0;
        data[isys].user_Av     = 30.0;
        data[isys].user_crflux = 2.5e-17;
        data[isys].user_GtoDN =
            (4.e0 * pi * rhoD * rD * rD * rD) / (3.e0 * DtoGM * amH);
    }

    Naunet naunet;
    naunet.Init();

    naunet.Reset(nsystem);

    double *y = new double[nsystem * NEQUATIONS];
    for (int isys = 0; isys < nsystem; isys++) {
        for (int i = 0; i < NEQUATIONS; i++) {
            y[isys * NEQUATIONS + i] = rawdata[isys][i + 2];
        }
    }

#ifdef NAUNET_DEBUG
    printf("Abundances in the first system\n");
    for (int i = 0; i < NEQUATIONS; i++) {
        printf("%13.7e ", y[i]);
    }
    printf("\n");
#endif

    FILE *fbin = fopen("evolution_paralleldata.bin", "w");
    FILE *ftxt = fopen("evolution_paralleldata.txt", "w");
    FILE *ttxt = fopen("time_paralleldata.txt", "w");

#ifdef NAUNET_DEBUG
    printf("Initialization is done. Start to evolve.\n");
    // FILE *rtxt = fopen("reactionrates.txt", "w");
    // double rates[NREACTIONS];
#endif

    double logtstart = 2.0, logtend = 5.0;
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

    naunet.Finalize();

    delete[] data;
    delete[] y;

    return 0;
}
