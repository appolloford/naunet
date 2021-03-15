#include <stdio.h>

#include "naunet.h"
#include "naunet_userdata.h"
#include "naunet_constants.h"
#include "fex.h"
#include "jtv.h"

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

    FILE *fbin = fopen("evolution.bin", "w");
    FILE *ftxt = fopen("evolution.txt", "w");
    realtype time = 0.0, dtyr = 1.0, tend = 1.e8;
    for (time = 0.0; time < tend; time += dtyr)
    {
        if (time < 1e5)
        {
            dtyr = fmax(9.0 * time, dtyr);
        }
        else
        {
            dtyr = 1e5;
        }
        naunet.solve(y, dtyr * spy, data);
        printf("Time = %13.7e yr\n", time);

        fwrite(&time, sizeof(realtype), 1, fbin);
        fwrite(y, sizeof(realtype), NSPECIES, fbin);

        fprintf(ftxt, "%13.7e ", time);
        for (int i = 0; i < NSPECIES; i++)
        {
            fprintf(ftxt, "%13.7e ", y[i]);
        }
        fprintf(ftxt, "\n");
    }

    fclose(fbin);
    fclose(ftxt);

    return 0;
}
