#include <stdio.h>

#include "naunet.h"
#include "naunet_userdata.h"
#include "naunet_macros.h"
#include "naunet_ode.h"
#include "naunet_timer.h"

int main()
{

    realtype spy = 86400.0 * 365.0;
    realtype pi = 3.14159265;
    realtype rD = 1.0e-5;
    realtype rhoD = 3.0;
    realtype DtoGM = 7.09e-3;
    realtype amH = 1.66043e-24;
    realtype nH = 1e5;
    realtype OPRH2 = 0.1;

    UserData *data = new UserData();
    data->nH = nH;
    data->Tgas = 15.0;
    data->user_Av = 30.0;
    data->user_crflux = 2.5e-17;
    data->user_GtoDN = (4.e0 * pi * rhoD * rD * rD * rD) / (3.e0 * DtoGM * amH);

    Naunet naunet;
    naunet.initSolver();

    realtype y[NSPECIES];
    for (int i = 0; i < NSPECIES; i++)
    {
        y[i] = 1.e-40;
    }
    y[IDX_pH2I] = 1.0 / (1.0 + OPRH2) * 0.5 * nH;
    y[IDX_oH2I] = OPRH2 / (1.0 + OPRH2) * 0.5 * nH;
    y[IDX_HDI] = 1.0e-5 * nH;
    y[IDX_HeI] = 1.0e-1 * nH;
    y[IDX_NI] = 2.1e-6 * nH;
    y[IDX_OI] = 1.8e-5 * nH;
    y[IDX_CI] = 7.3e-6 * nH;
    y[IDX_GRAIN0I] = 1.3215e-12 * nH;

    realtype time[10046];
    FILE *tfile = fopen("timeres.dat", "r");
    for (int i = 0; i < 10046; i++)
    {
        fscanf(tfile, "%lf\n", time + i);
    }
    fclose(tfile);

    FILE *fbin = fopen("evolution.bin", "w");
    FILE *ftxt = fopen("evolution.txt", "w");
#ifdef NAUNET_DEBUG
    FILE *rtxt = fopen("reactionrates.txt", "w");
    realtype rates[NREACTIONS];
#endif

    realtype dtyr = 1.0, tend = 1.e8;
    // for (time = 0.0; time < tend; time += dtyr)
    // {
    //     if (time < 1e5)
    //     {
    //         dtyr = fmax(9.0 * time, dtyr);
    //     }
    //     else
    //     {
    //         dtyr = 1e5;
    //     }
    for (int i = 0; i < 10046; i++)
    {

#ifdef NAUNET_DEBUG
        calculate_rates(rates, y, data);
        for (int j = 0; j < NREACTIONS; j++)
        {
            fprintf(rtxt, "%13.7e ", rates[j]);
        }
        fprintf(rtxt, "\n");
#endif

        dtyr = time[i + 1] - time[i];

        fwrite(time + i, sizeof(realtype), 1, fbin);
        fwrite(y, sizeof(realtype), NSPECIES, fbin);

        fprintf(ftxt, "%13.7e ", time[i]);
        for (int j = 0; j < NSPECIES; j++)
        {
            fprintf(ftxt, "%13.7e ", y[j]);
        }
        fprintf(ftxt, "\n");

        Timer timer;
        timer.start();
        naunet.solve(y, dtyr * spy, data);
        timer.stop();
        float duration = (float)timer.elapsed() / 1e6;
        printf("Time = %13.7e yr, elapsed: %8.5e sec\n", time[i + 1], duration);
    }

    fclose(fbin);
    fclose(ftxt);
#ifdef NAUNET_DEBUG
    fclose(rtxt);
#endif

    return 0;
}
